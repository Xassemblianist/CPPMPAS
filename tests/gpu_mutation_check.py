#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Show that the GPU tests (label gpu, tests/CMakeLists.txt) can fail: apply one one-line bug at a time to the source
of the CUDA build, rebuild, run the GPU tests the mutant must break, require each of them to fail (a test that did not
run, was skipped, timed out or passed does not count), restore the source. The tests must pass on the unmutated
source first, and the tests that killed a mutant must pass again on the restored source. Run it on the machine with
the GPU, in a CUDA build (-DCPPMPAS_ENABLE_CUDA=ON) whose CPU fixtures (initdump_*, dycore_*) the vs_cpu tests use.

usage: gpu_mutation_check.py <source dir> <build dir> [i,j,...]
The mutants: GPU-only bugs (the launch wrappers, kernel_fail, the stack limit, the skip path, device math, FMA
contraction, the CUDA flag policy) and mutants of tests/mutation_check.py, found by their fragment, for the GPU tests
of the dycore paths the GPU-only bugs do not single out (JW's Smagorinsky viscosity and 4th-order filters, the
restart run).
"""
import importlib.util
import os
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("mutation_check", HERE / "mutation_check.py")
mc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mc)

B, BC, X, L = "src/core/backend.hpp", "src/core/backend_cuda.cu", "tools/cppmpas.cpp", "src/core/libm.hpp"
CM = "CMakeLists.txt"
AS, DP = "src/core_atmosphere/atm_advance_scalars.cpp", "src/core_atmosphere/atm_dissipation_perturbation.cpp"
MP, KS = "src/core_atmosphere/atm_microphysics.cpp", "src/physics/mp_kessler.hpp"
RAD, RDRV = "src/physics/module_ra_rrtmg.cpp", "src/physics/mpas_atmphys_radiation.cpp"
RG_RAD = ["gpu_radiation_RG_0", "gpu_radiation_RG_1800"]


def dev(orig, mut):
    """the mutated statement in the device code only (the CPU build keeps the original): the vs_cpu tests see it.
    The mutated statement must use the same variables as the original: an extended lambda's implicit captures are
    the host pass's, so a variable only the device branch names is not captured (a first version that divided by
    rho_zz_old in place of rho_zz_new, or dropped pressure_b, read another capture's bytes on the device). It must
    also name them in the original's order: the device pass lays out the captures in the order of their first use
    in its body, so a swap of two captures reads each from the other's slot and cancels (a first o3clim mutant)."""
    return f"#if defined(__CUDA_ARCH__)\n{mut}\n#else\n{orig}\n#endif"


def SC(steps, *fields):
    return [f"gpu_dycore_SC_step{s}_{f}" for s in steps for f in fields]


SCV = ["gpu_dycore_SC_vs_cpu", "gpu_dycore_SC_mpas_compare"]
MW10 = [f"gpu_dycore_MW_step{s}_{f}" for s in (1, 10) for f in ("u", "theta_m")] + ["gpu_dycore_MW_vs_cpu"]
JWSTEP = [f"gpu_dycore_JW_step{s}_{f}" for s in (1, 10) for f in ("u", "theta_m")]
JW10 = JWSTEP + ["gpu_dycore_JW_vs_cpu"]
# (file, original fragment, mutated fragment, GPU tests that must fail)
GPU_MUTATIONS = [
    # a kernel skips the last element of its range
    (B, "  if (i <= long(last)) k(int(i));", "  if (i < long(last)) k(int(i));", ["gpu_device", "gpu_backend"]),
    # a whole-array kernel skips the last element (the garbage element of an MPAS array: only gpu_backend sees it)
    (B, "  if (i < n) k(i);", "  if (i + 1 < n) k(i);", ["gpu_backend"]),
    # the one-thread kernel does nothing
    (B, "template <class K> __global__ void kserial(K k) { k(); }",
     "template <class K> __global__ void kserial(K k) { (void)k; }", ["gpu_backend"] + MW10),
    # every launch reports a kernel failure
    (BC, "  if (g_err && g_err->set) {", "  if (g_err) {", ["gpu_device", "gpu_backend"]),
    # no device: exit 0 instead of the skip status
    (X, "      return gpu::kNoDeviceExit;", "      return 0;", ["gpu_skip_without_device"]),
    # device math: CUDA's pow in place of the libmvec lane on the device only (the host keeps glibc's)
    (L, "  return ::cppwrf::glibc::mvec::pow(arg(x), arg(y));",
     "#if defined(__CUDA_ARCH__)\n  return ::pow(x, y);\n#else\n  return ::cppwrf::glibc::mvec::pow(arg(x), arg(y));\n#endif",
     ["gpu_libm", "gpu_initdump_MW_exner", "gpu_initdump_JW_exner"] + MW10),
    # FMA contraction in the device code: nvcc's default --fmad=true (the strict flag left out)
    (CM, "$<$<COMPILE_LANGUAGE:CUDA>:--fmad=false;--ftz=false;", "$<$<COMPILE_LANGUAGE:CUDA>:--ftz=false;",
     ["cmake_cuda_strict_flags_last"] + MW10 + JW10),
    # kernel_fail does not end the thread: the statements after it run
    (B, '  detail::record_kernel_error(head, true, value, tail);\n  end_thread();',
     '  detail::record_kernel_error(head, true, value, tail);', ["gpu_backend"]),
    # a launch before anything initialised the runtime has no kernel error slot
    (B, "template <class Kern> void reserve_once(Kern* kernel) {\n  gpu::ensure_ready();",
     "template <class Kern> void reserve_once(Kern* kernel) {", ["gpu_backend"]),
    # the stack limit set to a kernel's frame even when it is smaller (lowered)
    (BC, "  if (a.localSizeBytes <= limit) {  // the limit is never lowered", "  if (false) {", ["gpu_backend"]),
    # the Release CUDA flags overwritten (the user's CMAKE_CUDA_FLAGS_RELEASE and CMake's -DNDEBUG dropped)
    (CM, "  enable_language(CUDA)\n", "  enable_language(CUDA)\n  set(CMAKE_CUDA_FLAGS_RELEASE \"-O3\")\n",
     ["cmake_cuda_strict_flags_last"]),
    # host fast-math through -Xcompiler not refused
    (CM, "|-Ofast|-f(fast-math|unsafe-math-optimizations|", "|-Ofast|-f(unsafe-math-optimizations|",
     ["cmake_rejects_cuda_6_CMAKE_CUDA_FLAGS", "cmake_rejects_cuda_6_CMAKE_CUDA_FLAGS_RELEASE"]),
    # --- the supercell's kernels, each mutant in the device code only ---
    # device math: CUDA's exp in place of the libmvec exp lane (Kessler's saturation pressure)
    (L, "CPPMPAS_FORCEINLINE double mvec_exp(double x) { return ::cppwrf::glibc::mvec::exp(arg(x)); }",
     "CPPMPAS_FORCEINLINE double mvec_exp(double x) {\n" + dev("  return ::cppwrf::glibc::mvec::exp(arg(x));",
                                                             "  return ::exp(x);") + "\n}", ["gpu_libm"]),
    # scalar transport, rk_steps 1-2 (TI:4772-4776): the 3rd-order weight's upwind sign flipped
    (AS, "        sw2(k, j) = m.adv_coefs(j, iEdge) + sign1(uhAvg(k, iEdge)) * m.adv_coefs_3rd(j, iEdge);",
     dev("        sw2(k, j) = m.adv_coefs(j, iEdge) + sign1(uhAvg(k, iEdge)) * m.adv_coefs_3rd(j, iEdge);",
         "        sw2(k, j) = m.adv_coefs(j, iEdge) - sign1(uhAvg(k, iEdge)) * m.adv_coefs_3rd(j, iEdge);"),
     SC((1, 10), "qv") + SCV),
    # monotonic transport, rk_step 3 (TI:5731-5735): the vertical flux divergence with the wrong sign of one face
    # (the divisor is rho_lim since the split transport, TI:5725-5735)
    (AS, "            (scalar_new(k, iCell) + (-(rdnw(k) * (wdtn(k + 1, iCell) - wdtn(k, iCell))))) / rho_lim(k, iCell);",
     dev("            (scalar_new(k, iCell) + (-(rdnw(k) * (wdtn(k + 1, iCell) - wdtn(k, iCell))))) / rho_lim(k, iCell);",
         "            (scalar_new(k, iCell) + (-(rdnw(k) * (wdtn(k + 1, iCell) + wdtn(k, iCell))))) / rho_lim(k, iCell);"),
     SC((1, 10), "qv") + SCV),
    # the limiter's bounds (TI:5374-5392) over five of the six neighbours: only the 30-min state sees it
    (AS, "        for (int j = 1; j <= 6; ++j) {\n          mx = vmax",
     dev("        for (int j = 1; j <= 6; ++j) {", "        for (int j = 1; j <= 5; ++j) {") + "\n          mx = vmax",
     ["gpu_dycore_SC_restart_qv", "gpu_dycore_SC_output30_qv"]),
    # perturbation mixing of u (DM:859-862): u + u_init in place of u - u_init
    (DP, "      u_mix[k] = u(k, iEdge) - u_init(k) * cosa + v_init(k) * sina;",
     dev("      u_mix[k] = u(k, iEdge) - u_init(k) * cosa + v_init(k) * sina;",
         "      u_mix[k] = u(k, iEdge) + u_init(k) * cosa + v_init(k) * sina;"), SC((10,), "u") + SCV),
    # perturbation mixing of theta (DM:1460-1462): t_init of the wrong level
    (DP, "              (((theta_m(k + 1, iCell) - t_init(k + 1, iCell)) - (theta_m(k, iCell) - t_init(k, iCell))) / (zp - z0) -",
     dev("              (((theta_m(k + 1, iCell) - t_init(k + 1, iCell)) - (theta_m(k, iCell) - t_init(k, iCell))) / (zp - z0) -",
         "              (((theta_m(k + 1, iCell) - t_init(k, iCell)) - (theta_m(k, iCell) - t_init(k, iCell))) / (zp - z0) -"),
     SC((10,), "theta_m") + SCV),
    # microphysics_to_MPAS (PI:903-934): theta_m from qc in place of qv
    (MP, "      theta_m(k, i) = P(th_p, i, k) * (1.0 + rvord * P(qv_p, i, k));",
     dev("      theta_m(k, i) = P(th_p, i, k) * (1.0 + rvord * P(qv_p, i, k));",
         "      theta_m(k, i) = P(th_p, i, k) * (1.0 + rvord * P(qc_p, i, k));"), SC((1, 10), "theta_m") + SCV),
    # the Kessler kernel on the device without libmvec's 2-lane functions (scalar pow / exp for every column)
    (KS, "CPPMPAS_HD inline bool lane(int i, int n) { return n > 1 && !((n & 1) && i == n - 1); }",
     "CPPMPAS_HD inline bool lane(int i, int n) {\n" + dev("  return n > 1 && !((n & 1) && i == n - 1);", "  return false;")
     + "\n}", ["gpu_dycore_SC_restart_qc", "gpu_dycore_SC_restart_qr", "gpu_dycore_SC_output30_qc",
               "gpu_dycore_SC_output30_qr"]),
    # the surface pressure of PI:950-964 with the base-state pressure of level 2
    (MP, "    surface_pressure(i) = surface_pressure(i) + pressure_p(1, i) + pressure_b(1, i);",
     dev("    surface_pressure(i) = surface_pressure(i) + pressure_p(1, i) + pressure_b(1, i);",
         "    surface_pressure(i) = surface_pressure(i) + pressure_p(1, i) + pressure_b(2, i);"),
     ["gpu_dycore_SC_restart_surface_pressure", "gpu_dycore_SC_output30_surface_pressure"]),
    # --- the radiation step of RG on the device (tests/radiation_replay), each mutant in the device code only ---
    # o3climatology_from_MPAS (OZ:317-321): the second month subtracted. (A first version swapped fact1 and fact2: it
    # survived, a mutant of nothing: the device pass orders an extended lambda's implicit captures by their first use
    # in the body it compiles, the host pass by its own, so a swap of two captures in the device branch reads each
    # from the other's slot and cancels. The mutated statement must use the captures in the original's order.)
    (RDRV, "    for (int k = 1; k <= noz; ++k) o3clim(k, iCell) = fact1 * ozmixm(nm, k, iCell) + fact2 * ozmixm(np, k, iCell);",
     dev("    for (int k = 1; k <= noz; ++k) o3clim(k, iCell) = fact1 * ozmixm(nm, k, iCell) + fact2 * ozmixm(np, k, iCell);",
         "    for (int k = 1; k <= noz; ++k) o3clim(k, iCell) = fact1 * ozmixm(nm, k, iCell) - fact2 * ozmixm(np, k, iCell);"),
     RG_RAD),
    # calc_cldfraction (CL:361): the exponent gamma in place of pexp. RG's 00:00 call has no cloud water or ice
    # anywhere (qcld < qcldmin in every cell, cldfrac = 0): only the 00:30 call reaches this statement
    (RDRV, "        cldfrac_p(i, k) = libm::pow(rhum / rhgrid, pexp) * (1. - libm::exp(arg));  // CL:361",
     dev("        cldfrac_p(i, k) = libm::pow(rhum / rhgrid, pexp) * (1. - libm::exp(arg));  // CL:361",
         "        cldfrac_p(i, k) = libm::pow(rhum / rhgrid, gamma) * (1. - libm::exp(arg));  // CL:361"), ["gpu_radiation_RG_1800"]),
    # the LW o3vmr diagnostic (RL:489): the pressures in units of 10 Pa
    (RDRV, "        p2dv(i, k) = pres_hyd_p(i, k) / 100.;",
     dev("        p2dv(i, k) = pres_hyd_p(i, k) / 100.;", "        p2dv(i, k) = pres_hyd_p(i, k) / 10.;"), RG_RAD),
    # rrtmg_lwrad's standard profile above the model top (LW 11993) matched at level kte + 1. (A first version asked for
    # at least two buffer layers, LW 11822: equivalent on RG, whose top pressures, 7.5-12.4 hPa, give 2-3 everywhere.)
    (RAD, "    v.tz[K(k)] = w.varint[K(k)] + (v.tz[K(kte)] - w.varint[K(kte)]);",
     dev("    v.tz[K(k)] = w.varint[K(k)] + (v.tz[K(kte)] - w.varint[K(kte)]);",
         "    v.tz[K(k)] = w.varint[K(k)] + (v.tz[K(kte + 1)] - w.varint[K(kte + 1)]);"), RG_RAD),
    # the LW stage of the level sums skips the top level of every column
    (RAD, "      if (lev > nlay[c]) return;",
     dev("      if (lev > nlay[c]) return;", "      if (lev >= nlay[c]) return;"), RG_RAD),
    # the SW zenith angle (SW 10206) with CUDA's cos in place of glibc's
    (RAD, "  const real coszrs = sinlat * sindec + coslat * cosdec * libm::cos(hrang);",
     dev("  const real coszrs = sinlat * sindec + coslat * cosdec * libm::cos(hrang);",
         "  const real coszrs = sinlat * sindec + coslat * cosdec * ::cos(hrang);"), RG_RAD),
    # the SW McICA subcolumns (SW 10540-10549) with permuteseed 2 (no cloud at 00:00: the 00:30 call)
    (RAD, "  mcica_subcol(icld, 1, 112, nlay, v.pavel, v.cldfrac, w.hlay, in.idcor, in.julday, lat,",
     dev("  mcica_subcol(icld, 1, 112, nlay, v.pavel, v.cldfrac, w.hlay, in.idcor, in.julday, lat,",
         "  mcica_subcol(icld, 2, 112, nlay, v.pavel, v.cldfrac, w.hlay, in.idcor, in.julday, lat,"), ["gpu_radiation_RG_1800"]),
]
# tests/mutation_check.py mutants by a piece of their original fragment, with the GPU tests they must break (a mutant
# of the shared source changes the CPU build alike: the vs_cpu tests are not among them)
CPU_MUTANTS = [
    # JW: smagorinsky_2d's eddy viscosity from u twice instead of u and v
    ("smagorinsky_2d(m, f.eddy_visc_horz, f.u2, d.v,", JWSTEP),
    # JW: the del4 vorticity of u_dissipation_3d (DM:735-751)
    ("        const real edge_sign = m.invAreaTriangle(iVertex) * m.dcEdge(iEdge) * m.edgesOnVertex_sign(i, iVertex);",
     JWSTEP),
    # the RBF reconstruction coefficients and the local vertical unit vectors (init block, not in the initdump set)
    ("    alpha = alpha / static_cast<real>(pointCount);", ["gpu_init_restart_MW", "gpu_init_restart_JW"]),
    ("      mpas_unit_vec_in_r3(v);\n      for (int d = 1; d <= 3; ++d) lvu(d, iCell) = v[d - 1];",
     ["gpu_init_restart_JW"]),
    # the restart run: the coupled fields of the restart stream recomputed instead of kept (mpas_atm_core.F:522)
    ("  if (!cfg.config_do_restart || (cfg.config_do_restart && cfg.config_do_DAcycling)) atm_in",
     ["gpu_dycore_restart0_MW_rho_zz", "gpu_dycore_restart0_MW_exner", "gpu_dycore_restart0_MW_ru"]),
]


def cpu_mutant(piece):
    hits = [m for m in mc.MUTATIONS if piece in m[1]]
    if len(hits) != 1:
        raise SystemExit(f"tests/mutation_check.py: {len(hits)} mutants with the fragment {piece!r}")
    return hits[0][:3]


# --- the physics on the GPU (RG): each mutant in the device code only, appended after the indices above ---
PH, DG = "src/physics/", "src/core_atmosphere/diagnostics/"
RGR = lambda *f: [f"gpu_dycore_RG_restart_{x}" for x in f]
RGV = ["gpu_dycore_RG_vs_cpu"]
RGS = lambda f: [f"gpu_dycore_RG_streams_{f}"]
RG24 = lambda f: [f"gpu_dycore_RG24_{f}"]
PHYSICS_MUTATIONS = [
    # MPAS_to_physics' model-top pressure (PI:494-506) without the libmvec log / exp lanes
    (PH + "mpas_atmphys_interface.cpp", "      if (lane(i - its, n))",
     dev("      if (lane(i - its, n))", "      if (false && lane(i - its, n))"), RGR("theta_m", "hfx") + RGV),
    # sfclayer_from_MPAS (SL:422-470): twice the roughness length
    (PH + "mpas_atmphys_sfclayer.cpp", "    p.znt[I(i)] = m.znt(i);",
     dev("    p.znt[I(i)] = m.znt(i);", "    p.znt[I(i)] = 2. * m.znt(i);"), RGR("ust", "hfx") + RGV),
    # Noah's SFLX (ICE = 0 points) skipped
    (PH + "module_sf_noahdrv.cpp", "    nk(ii, 0);   // returns at once unless ICE = 0 on land (active)",
     dev("    nk(ii, 0);   // returns at once unless ICE = 0 on land (active)", "    (void)nk;"),
     RGR("skintemp", "tslb") + RGV),
    # the sea-ice driver's SFLX skipped
    (PH + "module_sf_noahdrv.cpp", "  parallel_for(a.its, a.ite, [=] CPPMPAS_LAMBDA(int i) { si(i - ims, 0); });",
     "  parallel_for(a.its, a.ite, [=] CPPMPAS_LAMBDA(int i) {\n" + dev("    si(i - ims, 0);", "    (void)si;")
     + "\n  });", RGR("skintemp", "hfx")),
    # seaice_to_MPAS (seaice.F:387-405): open-water albedo 0.07
    (PH + "mpas_atmphys_driver_lsm.cpp", "      L.at(L.sfc_albedo_p, i) = x * L.at(L.sfc_albedo_p, i) + (1. - x) * 0.08;",
     dev("      L.at(L.sfc_albedo_p, i) = x * L.at(L.sfc_albedo_p, i) + (1. - x) * 0.08;",
         "      L.at(L.sfc_albedo_p, i) = x * L.at(L.sfc_albedo_p, i) + (1. - x) * 0.07;"), RGR("sfc_albedo")),
    # atmphys_sfc_diagnostics: th2m one kelvin off
    (PH + "mpas_atmphys_driver_lsm.cpp", "    p.th2m(i) = p.t2m(i) * libm::pow(pc::P0 / p.surface_pressure(i), pc::rcp);",
     dev("    p.th2m(i) = p.t2m(i) * libm::pow(pc::P0 / p.surface_pressure(i), pc::rcp);",
         "    p.th2m(i) = p.t2m(i) * libm::pow(pc::P0 / p.surface_pressure(i), pc::rcp) + 1.;"), RGR("th2m")),
    # YSU's column kernel skipped
    (PH + "mpas_atmphys_pbl.cpp", "    parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) { k(i - 1, 0); });",
     "    parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {\n" + dev("      k(i - 1, 0);", "      (void)k;") + "\n    });",
     RGR("rublten", "rthblten") + RGV),
    # GWDO's column routine skipped
    (PH + "mpas_atmphys_gwdo.cpp", "    bl_gwdo::bl_gwdo_run<real>(c, kte, gravity, cp, R_d, ep_1, pii, m.dt_pbl);",
     dev("    bl_gwdo::bl_gwdo_run<real>(c, kte, gravity, cp, R_d, ep_1, pii, m.dt_pbl);", "    (void)c;"),
     RGR("dtaux3d", "dusfcg")),
    # New Tiedtke's cuinin without the libmvec lanes
    (PH + "module_cu_ntiedtke.cpp", "    c.cuinin_mvec = lane(i - its, im);",
     dev("    c.cuinin_mvec = lane(i - its, im);", "    c.cuinin_mvec = false && lane(i - its, im);"),
     RGR("rthcuten", "rqvcuten")),
    # CV:843-849: the forcing's temperature tendency without the exner factor
    (PH + "mpas_atmphys_convection.cpp",
     "          (m.rthdynten(k, i) + m.rthblten(k, i) + m.rthratenlw(k, i) + m.rthratensw(k, i)) * v.pi_p(i, k);",
     dev("          (m.rthdynten(k, i) + m.rthblten(k, i) + m.rthratenlw(k, i) + m.rthratensw(k, i)) * v.pi_p(i, k);",
         "          (m.rthdynten(k, i) + m.rthblten(k, i) + m.rthratenlw(k, i) + m.rthratensw(k, i)) * 1.0;"),
     RGR("rthcuten", "rqvcuten")),
    # WSM6 without the libmvec lanes of the warm-rain loop (every column the scalar log / exp)
    (PH + "module_mp_wsm6.cpp", "           sr, i - its, ite - its + 1);",
     dev("           sr, i - its, ite - its + 1);", "           sr, 0 * (i - its), 1 + 0 * (ite - its + 1));"),
     RGR("qr", "qv") + RGV),
    # compute_radar_reflectivity skipped (l_diags)
    (PH + "mpas_atmphys_microphysics.cpp",
     "      wsm6::refl10cm_wsm6(*m.radar, qv1d, qr1d, qs1d, qg1d, t1d, p1d, dBZ1d, kte);",
     dev("      wsm6::refl10cm_wsm6(*m.radar, qv1d, qr1d, qs1d, qg1d, t1d, p1d, dBZ1d, kte);", "      (void)m.radar;"),
     RG24("diag.2026-10-12_03.00.00")),
    # compute_relhum (DRV:798-850) in tenths
    (PH + "mpas_atmphys_microphysics.cpp", "        m.relhum(k, i) = m.relhum(k, i) * 100.;",
     dev("        m.relhum(k, i) = m.relhum(k, i) * 100.;", "        m.relhum(k, i) = m.relhum(k, i) * 10.;"),
     RGS("history.2026-10-12_01.00.00")),
    # physics_get_tend (TD:382-426): the PBL's edge tendency subtracted
    (PH + "mpas_atmphys_todynamics.cpp",
     "      for (int k = 1; k <= nV; ++k) tend_u(k, i) = tend_u(k, i) + m.rublten_Edge(k, i) * mass_edge(k, i);",
     dev("      for (int k = 1; k <= nV; ++k) tend_u(k, i) = tend_u(k, i) + m.rublten_Edge(k, i) * mass_edge(k, i);",
         "      for (int k = 1; k <= nV; ++k) tend_u(k, i) = tend_u(k, i) - m.rublten_Edge(k, i) * mass_edge(k, i);"),
     ["gpu_dycore_RG_step1_u"] + RGV),
    # TI:2777-2783: rthdynten without the qv term
    (PH + "mpas_atmphys_driver.cpp",
     "      rthdynten(k, iCell) = fac_m * (rthdynten(k, iCell) - theta_local * rv / rgas * rqvdynten(k, iCell));",
     dev("      rthdynten(k, iCell) = fac_m * (rthdynten(k, iCell) - theta_local * rv / rgas * rqvdynten(k, iCell));",
         "      rthdynten(k, iCell) = fac_m * (rthdynten(k, iCell) + 0. * theta_local * rv / rgas * rqvdynten(k, iCell));"),
     RGR("rthdynten")),
    # interp_tofixed_pressure: the trapping level one too high
    (DG + "mpas_isobaric_diagnostics.cpp", "          kupper[icol] = kk;",
     dev("          kupper[icol] = kk;", "          kupper[icol] = kk > 1 ? kk - 1 : kk;"),
     RGS("diag.2026-10-12_00.00.00")),
    # the updraft helicity's vertex-to-cell sums without the kite areas
    (DG + "mpas_convective_diagnostics.cpp",
     "            updraft_helicity(k, iCell) = updraft_helicity(k, iCell) + kiteAreasOnVertex(j, i) * vorticity(k, i);",
     dev("            updraft_helicity(k, iCell) = updraft_helicity(k, iCell) + kiteAreasOnVertex(j, i) * vorticity(k, i);",
         "            updraft_helicity(k, iCell) = updraft_helicity(k, iCell) + 0. * kiteAreasOnVertex(j, i) * vorticity(k, i);"),
     RG24("diag.2026-10-12_03.00.00")),
    # getcape skipped
    (DG + "mpas_convective_diagnostics.cpp",
     "      getcape(nVertLevels, p_in.data(), t_in.data(), td_in.data(), cape_out, cin_out);",
     dev("      getcape(nVertLevels, p_in.data(), t_in.data(), td_in.data(), cape_out, cin_out);", "      (void)p_in;"),
     RGS("diag.2026-10-12_00.00.00")),
    # KernelCopy in the CUDA build copies zero tables in place of the scheme's (the CPU build reads the tables directly)
    (PH + "mpas_atmphys_interface.hpp", "  void set(const T& x) { s.assign(1, x); }", "  void set(const T& x) { (void)x; s.assign(1, T{}); }",
     ["gpu_dycore_RG_step1_theta_m"] + RGV),
]


def mutants():
    return [m for m in GPU_MUTATIONS] + [cpu_mutant(p) + (t,) for p, t in CPU_MUTANTS] + PHYSICS_MUTATIONS


def build(b):
    # -j4 on the shared workstation; GPU_MUT_JOBS on a machine of its own
    return subprocess.run(["cmake", "--build", str(b), "-j" + os.environ.get("GPU_MUT_JOBS", "4")],
                          capture_output=True, text=True)


def main():
    if len(sys.argv) not in (3, 4):
        raise SystemExit(__doc__)
    src, b = Path(sys.argv[1]), Path(sys.argv[2])
    muts = mutants()
    if len(sys.argv) == 4:
        muts = [muts[int(i)] for i in sys.argv[3].split(",")]
    bad, covered = [], set()
    r0 = build(b)
    if r0.returncode != 0:
        raise SystemExit("the unmutated source does not build")
    needed = sorted({t for *_, ts in muts for t in ts})
    r = mc.ctest(b, mc.any_of(needed))
    broken = mc.not_passing(r.stdout, needed)
    for t in sorted(broken):
        bad.append(f"{t} does not pass on the unmutated source (status {mc.test_status(r.stdout, t)!r})")
    for path, old, new, tests in muts:
        p = src / path
        orig = p.read_text()
        if orig.count(old) != 1:
            bad.append(f"{path}: fragment found {orig.count(old)} times: {old!r}")
            continue
        try:
            p.write_text(orig.replace(old, new))
            rb = build(b)
            if rb.returncode != 0:
                bad.append(f"{path}: mutant does not build: {new.strip()[:80]}")
                print(rb.stdout[-1500:], rb.stderr[-1500:])
                continue
            r = mc.ctest(b, mc.any_of(tests))
            for t in tests:
                status = mc.test_status(r.stdout, t)
                k = mc.is_kill(r.returncode, status, t, broken)
                print(f"{'KILLED' if k else 'SURVIVED'} {t:40s} by {path}: {new.strip()[:70]!r}", flush=True)
                if k:
                    covered.add(t)
                else:
                    bad.append(f"{t} survived {path}: {new.strip()[:80]} (status {status!r})")
        finally:
            p.write_text(orig)
    if build(b).returncode != 0:
        bad.append("the restored source does not build")
    if covered:
        r = mc.ctest(b, mc.any_of(sorted(covered)))
        for t in sorted(mc.not_passing(r.stdout, covered)):
            bad.append(f"{t} does not pass on the restored source: its kills are not trusted")
    for x in bad:
        print("PROBLEM", x)
    print(f"{len(covered)} GPU tests killed by at least one mutant; {len(bad)} problems")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
