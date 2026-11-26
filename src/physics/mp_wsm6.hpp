// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// WSM6 cloud microphysics (MMM-physics mp_wsm6.F90, pinned commit a4baf7f), the column scheme of CPPWRF (commit
// 7fac5bf, src/phys/mp_wsm6/wsm6.hpp lines 16-2106, last changed in 802ce8c): `Wsm6Params` (mp_wsm6_init and,
// in operator(), mp_wsm6_run for one column). Shared with CPPWRF (CONSTITUTION.md 2.6, SHARED_WITH_CPPWRF.md).
//
// Same code: MPAS v8.4.2 builds WSM6 from MMM-physics a4baf7f (the cppmpas-data `inputs` asset
// mmm-physics-a4baf7f.tar.gz, sha256 a7bbe605...); its mp_wsm6.F90, mp_wsm6_effectRad.F90, mp_radar.F90 and
// module_libmassv.F90 are byte-identical to WRF 4.8.0's phys/physics_mmm files that CPPWRF ported (diff: none).
// MPAS's wrapper physics_wrf/module_mp_wsm6.F differs from WRF 4.8.0's phys/module_mp_wsm6.F only in kinds
// (RKIND for kind_phys) and in not calling refl10cm_wsm6 (MPAS computes the reflectivity in the driver at
// diagnostic times, mpas_atmphys_driver_microphysics.F:637).
//
// The struct is CPPWRF's, copied unchanged except
//  * CPPWRF_HD spelled CPPMPAS_HD; libm:: resolves to cppmpas::libm (src/core/libm.hpp: glibc's __exp_fma,
//    __log_fma, __ieee754_pow_fma, __atan_fma, as CPPWRF's cppwrf::libm);
//  * the MPAS parameter lane_i / lane_n of operator() (default 0 / 1: CPPWRF's kernel). MPAS builds mp_wsm6.F90
//    with -O3 where WRF uses -O2 -ftree-vectorize, and gfortran -O3 vectorises one more loop: work1 / work2,
//    mp_wsm6.F90:850-856 (diffac, diffac, venfac). Its 4 log and 4 exp per iteration become libmvec 2-lane calls
//    (_ZGVbN2v_log / _ZGVbN2v_exp, lines 852-854), in pairs of i from its; an odd last i, and a 1-trip loop, use
//    the scalar log@plt / exp@plt. MPAS calls operator() for column i of its..ite with lane_i = i - its,
//    lane_n = ite - its + 1 (vm::lane, CPPWRF's lane(i, n) rule); the warm-rain loop computes these values as
//    diffac_v / venfac_v. Every other libm call of mp_wsm6_run, slope_wsm6, nislfv_rain_plm(6) and mp_wsm6_init
//    is scalar in the MPAS executable as in WRF's.
//    Object code: mp_wsm6.F90 compiled standalone with the reference flags (gfortran 14.2.0 -O3
//    -ffree-line-length-none -fconvert=big-endian -ffree-form -fdefault-real-8 -fdefault-double-8) is
//    instruction-identical to the reference executable's (bin/atmosphere_model, sha256 524a1cf4...,
//    __mp_wsm6_MOD_mp_wsm6_run 0x499780, slope_wsm6 0x493780, nislfv_rain_plm 0x493fd0, nislfv_rain_plm6
//    0x495b40, mp_wsm6_init 0x4a9ce0; addresses and padding aside); with -g it maps the four _ZGVbN2v_log and four
//    _ZGVbN2v_exp calls of mp_wsm6_run to lines 852-854 and the scalar copies to the same lines. The vector loop
//    is entered behind a 1-trip test and overlap tests of work1 / work2 with t, p, den, which the wrapper's local
//    copies never meet. No FMA, no ymm register in any of these functions.
//
// CPPMPAS is FP64 only (core/precision.hpp): the float paths and the CPPWRF_DOUBLE stage helpers of the copy are
// not compiled. tests/wsm6_column checks the kernel against MPAS's own WSM6 built with the reference flags.
#include <algorithm>
#include <cassert>
#include <cmath>
#include <type_traits>

#include "core/backend.hpp"
#include "core/libm.hpp"

namespace cppmpas::phys::wsm6 {

// CPPWRF's lane rule (CPPWRF 7fac5bf src/core/backend.hpp:197): libmvec's 2-lane function for the pairs of
// iterations from the start of the loop, the scalar libm for an odd last iteration and a 1-trip loop
namespace vm {
CPPMPAS_HD inline bool lane(int i, int n) { return n > 1 && !((n & 1) && i == n - 1); }
CPPMPAS_HD inline double exp(double a, int i, int n) {
  return lane(i, n) ? ::cppmpas::libm::mvec_exp(a) : ::cppmpas::libm::exp(a);
}
CPPMPAS_HD inline double log(double a, int i, int n) {
  return lane(i, n) ? ::cppmpas::libm::mvec_log(a) : ::cppmpas::libm::log(a);
}
}  // namespace vm

// ---- copied from CPPWRF 7fac5bf src/phys/mp_wsm6/wsm6.hpp:16-2106 (MPAS parameter marked "MPAS parameter") ----
constexpr int kMaxLevels = 128;

// The type of the floating literals of the per-step code: Real (float in the FP32 builds: WRF real*4's default-kind
// literals, CPU and CUDA alike; an unsuffixed C++ literal is double and would turn the arithmetic into FP64; double
// in the FP64 builds).
template <typename Real>
using wlit = Real;

// rgmma (mp_wsm6.F90, module_mp_wsm5_accel.F, module_mp_wdm5.F / wdm6): real(kind_phys) / REAL, single precision
// in WRF real*4 (ccpp_kind_types: kind_phys = selected_real_kind(6) without DOUBLE_PRECISION), the product of the
// 10000 factors in float; the FP64 builds below as before
template <typename Real>
inline Real rgmma_calc(Real x) {
    if (x == static_cast<Real>(1.0)) return static_cast<Real>(0.0);
    if constexpr (std::is_same_v<Real, float>) {
        const float euler = 0.577215664901532f;  // a real(kind_phys) parameter
        float r = x * libm::exp(euler * x);
        for (int i = 1; i <= 10000; ++i) {
            const float y = float(i);
            r = r * (1.000f + x / y) * libm::exp(-x / y);
        }
        return 1.f / r;
    }
    constexpr double euler = 0.577215664901532;
    double r = static_cast<double>(x) * libm::exp(euler * static_cast<double>(x));
    for (int i = 1; i <= 10000; ++i) {
        double y = static_cast<double>(i);
        r = r * (1.000 + static_cast<double>(x) / y) * libm::exp(-static_cast<double>(x) / y);
    }
    return static_cast<Real>(1.0 / r);
}

template <typename Real>
struct Wsm6Params {
    Real g{9.81};
    Real cpd{1004.5};
    Real cpv{1846.4};
    Real rd{287.0};
    Real rv{461.6};
    Real t0c{273.15};
    Real ep1{0.6083623693379791};
    Real ep2{0.6217504332755632};
    Real qmin{1.0e-15};
    Real xls{2.85e6};
    Real xlv0{2.50e6};
    Real xlf0{0.35e6};
    Real den0{1.28};
    Real denr{1000.0};
    Real dens{100.0};
    Real cliq{4190.0};
    Real cice{2106.0};
    Real psat{610.78};
    int hail_opt{0};

    Real dtcldcr{120.0};
    Real n0r{8.0e6};
    Real avtr{841.9};
    Real bvtr{0.8};
    Real r0{0.8e-5};
    Real peaut{0.55};
    Real xncr{3.0e8};
    Real xmyu{1.718e-5};
    Real avts{11.72};
    Real bvts{0.41};
    Real lamdarmax{8.0e4};
    Real lamdasmax{1.0e5};
    Real dicon{11.9};
    Real dimax{500.0e-6};
    Real pfrz1{100.0};
    Real pfrz2{0.66};
    Real qcrmin{1.0e-9};
    Real eacrc{1.0};
    Real qs0{6.0e-4};
    Real n0smax{1.0e11};
    Real n0s{2.0e6};
    Real alpha{0.12};

    Real n0g{4.0e6};
    Real deng{500.0};
    Real avtg{330.0};
    Real bvtg{0.8};
    Real lamdagmax{6.0e4};

    Real pi{3.14159265358979323846};
    Real xlv1{0.0};
    Real qc0{0.0};
    Real qck1{0.0};
    Real pidnc{0.0};

    Real bvtr1{0.0}, bvtr2{0.0}, bvtr3{0.0}, bvtr4{0.0}, bvtr6{0.0};
    Real g1pbr{0.0}, g3pbr{0.0}, g4pbr{0.0}, g6pbr{0.0}, g5pbro2{0.0};
    Real pvtr{0.0}, eacrr{1.0}, pacrr{0.0}, precr1{0.0}, precr2{0.0}, roqimax{0.0};

    Real bvts1{0.0}, bvts2{0.0}, bvts3{0.0}, bvts4{0.0};
    Real g1pbs{0.0}, g3pbs{0.0}, g4pbs{0.0}, g5pbso2{0.0};
    Real pvts{0.0}, pacrs{0.0}, precs1{0.0}, precs2{0.0}, pidn0r{0.0}, pidn0s{0.0}, pacrc{0.0};

    Real bvtg1{0.0}, bvtg2{0.0}, bvtg3{0.0}, bvtg4{0.0};
    Real g1pbg{0.0}, g3pbg{0.0}, g4pbg{0.0}, g5pbgo2{0.0};
    Real pvtg{0.0}, pacrg{0.0}, precg1{0.0}, precg2{0.0}, pidn0g{0.0};

    Real rslopermax{0.0}, rslopesmax{0.0}, rslopegmax{0.0};
    Real rsloperbmax{0.0}, rslopesbmax{0.0}, rslopegbmax{0.0};
    Real rsloper2max{0.0}, rslopes2max{0.0}, rslopeg2max{0.0};
    Real rsloper3max{0.0}, rslopes3max{0.0}, rslopeg3max{0.0};

    void init(Real p_den0, Real p_denr, Real p_dens, Real p_cliq, Real p_cpv, int p_hail_opt) {
        den0 = p_den0;
        denr = p_denr;
        dens = p_dens;
        cliq = p_cliq;
        cpv = p_cpv;
        hail_opt = p_hail_opt;

        if (hail_opt == 1) {
            n0g = 4.0e4;
            deng = 700.0;
            avtg = 285.0;
            bvtg = 0.8;
            lamdagmax = 2.0e4;
        } else {
            n0g = 4.0e6;
            deng = 500.0;
            avtg = 330.0;
            bvtg = 0.8;
            lamdagmax = 6.0e4;
        }

        pi = 4.0 * libm::atan(1.0);
        xlv1 = cliq - cpv;

        qc0 = (4.0 / 3.0) * pi * denr * (r0 * r0 * r0) * xncr / den0;
        qck1 = 0.104 * 9.8 * peaut / libm::pow(xncr * denr, 1.0 / 3.0) / xmyu * libm::pow(den0, 4.0 / 3.0);
        pidnc = pi * denr / 6.0;

        bvtr1 = 1.0 + bvtr;
        bvtr2 = 2.5 + 0.5 * bvtr;
        bvtr3 = 3.0 + bvtr;
        bvtr4 = 4.0 + bvtr;
        bvtr6 = 6.0 + bvtr;
        g1pbr = rgmma_calc(bvtr1);
        g3pbr = rgmma_calc(bvtr3);
        g4pbr = rgmma_calc(bvtr4);
        g6pbr = rgmma_calc(bvtr6);
        g5pbro2 = rgmma_calc(bvtr2);
        pvtr = avtr * g4pbr / 6.0;
        eacrr = 1.0;
        pacrr = pi * n0r * avtr * g3pbr * 0.25 * eacrr;
        precr1 = 2.0 * pi * n0r * 0.78;
        precr2 = 2.0 * pi * n0r * 0.31 * std::sqrt(avtr) * g5pbro2;
        roqimax = 2.08e22 * libm::pow(dimax, 8.0);

        bvts1 = 1.0 + bvts;
        bvts2 = 2.5 + 0.5 * bvts;
        bvts3 = 3.0 + bvts;
        bvts4 = 4.0 + bvts;
        g1pbs = rgmma_calc(bvts1);
        g3pbs = rgmma_calc(bvts3);
        g4pbs = rgmma_calc(bvts4);
        g5pbso2 = rgmma_calc(bvts2);
        pvts = avts * g4pbs / 6.0;
        pacrs = pi * n0s * avts * g3pbs * 0.25;
        precs1 = 4.0 * n0s * 0.65;
        precs2 = 4.0 * n0s * 0.44 * std::sqrt(avts) * g5pbso2;
        pidn0r = pi * denr * n0r;
        pidn0s = pi * dens * n0s;
        pacrc = pi * n0s * avts * g3pbs * 0.25 * eacrc;

        bvtg1 = 1.0 + bvtg;
        bvtg2 = 2.5 + 0.5 * bvtg;
        bvtg3 = 3.0 + bvtg;
        bvtg4 = 4.0 + bvtg;
        g1pbg = rgmma_calc(bvtg1);
        g3pbg = rgmma_calc(bvtg3);
        g4pbg = rgmma_calc(bvtg4);
        g5pbgo2 = rgmma_calc(bvtg2);
        pvtg = avtg * g4pbg / 6.0;
        pacrg = pi * n0g * avtg * g3pbg * 0.25;
        precg1 = 2.0 * pi * n0g * 0.78;
        precg2 = 2.0 * pi * n0g * 0.31 * std::sqrt(avtg) * g5pbgo2;
        pidn0g = pi * deng * n0g;

        rslopermax = 1.0 / lamdarmax;
        rslopesmax = 1.0 / lamdasmax;
        rslopegmax = 1.0 / lamdagmax;
        rsloperbmax = libm::pow(rslopermax, bvtr);
        rslopesbmax = libm::pow(rslopesmax, bvts);
        rslopegbmax = libm::pow(rslopegmax, bvtg);
        rsloper2max = rslopermax * rslopermax;
        rslopes2max = rslopesmax * rslopesmax;
        rslopeg2max = rslopegmax * rslopegmax;
        rsloper3max = rsloper2max * rslopermax;
        rslopes3max = rslopes2max * rslopesmax;
        rslopeg3max = rslopeg2max * rslopegmax;
        if constexpr (std::is_same_v<Real, float>) {
            // mp_wsm6_init in WRF real*4: the products above in REAL from left to right (the double evaluation
            // rounds qck1, roqimax and precg2, with hail_opt 1 precg1, to other floats); products of literals and
            // parameters folded in REAL, parameter**.5 and dimax**8 folded correctly rounded, avtg**.5 a runtime powf
            const float sqavtr = std::sqrt(avtr), sqavts = std::sqrt(avts);
            qc0 = 4.f / 3.f * pi * denr * (r0 * r0 * r0) * xncr / den0;
            qck1 = 0.104f * 9.8f * peaut / libm::pow(xncr * denr, 1.f / 3.f) / xmyu * libm::pow(den0, 4.f / 3.f);
            pidnc = pi * denr / 6.f;
            pvtr = avtr * g4pbr / 6.f;
            pacrr = pi * n0r * avtr * g3pbr * .25f * eacrr;
            precr1 = 2.f * pi * n0r * .78f;
            precr2 = 2.f * pi * n0r * .31f * sqavtr * g5pbro2;
            roqimax = 2.08e22f * static_cast<float>(libm::pow(static_cast<double>(dimax), 8.0));
            pvts = avts * g4pbs / 6.f;
            pacrs = pi * n0s * avts * g3pbs * .25f;
            precs1 = 4.f * n0s * .65f;
            precs2 = 4.f * n0s * .44f * sqavts * g5pbso2;
            pidn0r = pi * denr * n0r;
            pidn0s = pi * dens * n0s;
            pacrc = pi * n0s * avts * g3pbs * .25f * eacrc;
            pacrg = pi * n0g * avtg * g3pbg * .25f;
            pvtg = avtg * g4pbg / 6.f;
            precg1 = 2.f * pi * n0g * .78f;
            precg2 = 2.f * pi * n0g * .31f * libm::pow(avtg, .5f) * g5pbgo2;
            pidn0g = pi * deng * n0g;
        }
    }

    CPPMPAS_HD inline Real cpmcal(Real x) const {
        return cpd * (wlit<Real>(1.0) - std::max(x, qmin)) + std::max(x, qmin) * cpv;
    }

    CPPMPAS_HD inline Real xlcal(Real x) const {
        return xlv0 - xlv1 * (x - t0c);
    }

    CPPMPAS_HD inline Real diffus(Real x, Real y) const {
        return wlit<Real>(8.794e-5) * libm::exp(libm::log(x) * wlit<Real>(1.81)) / y;
    }

    CPPMPAS_HD inline Real viscos(Real x, Real y) const {
        return wlit<Real>(1.496e-6) * (x * std::sqrt(x)) / (x + wlit<Real>(120.0)) / y;
    }

    CPPMPAS_HD inline Real xka(Real x, Real y) const {
        return wlit<Real>(1.414e3) * viscos(x, y) * y;
    }

    CPPMPAS_HD inline Real diffac(Real a, Real b, Real c, Real d, Real e) const {
        return d * a * a / (xka(c, d) * rv * c * c) + wlit<Real>(1.0) / (e * diffus(c, b));
    }

    CPPMPAS_HD inline Real venfac(Real a, Real b, Real c) const {
        Real v = viscos(b, c);
        return libm::exp(libm::log(v / diffus(b, a)) * wlit<Real>(0.3333333)) / std::sqrt(v) * std::sqrt(std::sqrt(den0 / c));
    }

    // MPAS parameter: diffus / diffac / venfac of the work1 / work2 loop (mp_wsm6.F90:850-856) with the lanes of
    // MPAS's -O3 build (vm::lane(lane_i, lane_n)); with lane_n = 1 they are diffac / venfac above
    CPPMPAS_HD inline Real diffus_v(Real x, Real y, int li, int ln) const {
        return wlit<Real>(8.794e-5) * vm::exp(vm::log(x, li, ln) * wlit<Real>(1.81), li, ln) / y;
    }
    CPPMPAS_HD inline Real diffac_v(Real a, Real b, Real c, Real d, Real e, int li, int ln) const {
        return d * a * a / (xka(c, d) * rv * c * c) + wlit<Real>(1.0) / (e * diffus_v(c, b, li, ln));
    }
    CPPMPAS_HD inline Real venfac_v(Real a, Real b, Real c, int li, int ln) const {
        Real v = viscos(b, c);
        return vm::exp(vm::log(v / diffus_v(b, a, li, ln), li, ln) * wlit<Real>(0.3333333), li, ln) / std::sqrt(v)
            * std::sqrt(std::sqrt(den0 / c));
    }

    CPPMPAS_HD inline Real conden(Real a, Real b, Real c, Real d, Real e) const {
        return (std::max(b, qmin) - c) / (wlit<Real>(1.0) + d * d / (rv * e) * c / (a * a));
    }

    CPPMPAS_HD inline void slope_wsm6(int km, const Real* qr, const Real* qs, const Real* qg,
                                    const Real* den, const Real* denfac, const Real* t,
                                    Real rslope[][3], Real rslopeb[][3], Real rslope2[][3],
                                    Real rslope3[][3], Real vt[][3]) const {
        for (int k = 0; k < km; ++k) {
            Real supcol = t0c - t[k];
            Real n0sfac = std::max(std::min(libm::exp(alpha * supcol), n0smax / n0s), static_cast<Real>(1.0));

            // Rain
            if (qr[k] <= qcrmin) {
                rslope[k][0] = rslopermax;
                rslopeb[k][0] = rsloperbmax;
                rslope2[k][0] = rsloper2max;
                rslope3[k][0] = rsloper3max;
            } else {
                Real lamdar = std::sqrt(std::sqrt(pidn0r / (qr[k] * den[k])));
                rslope[k][0] = wlit<Real>(1.0) / lamdar;
                rslopeb[k][0] = libm::pow(rslope[k][0], bvtr);
                rslope2[k][0] = rslope[k][0] * rslope[k][0];
                rslope3[k][0] = rslope2[k][0] * rslope[k][0];
            }

            // Snow
            if (qs[k] <= qcrmin) {
                rslope[k][1] = rslopesmax;
                rslopeb[k][1] = rslopesbmax;
                rslope2[k][1] = rslopes2max;
                rslope3[k][1] = rslopes3max;
            } else {
                Real lamdas = std::sqrt(std::sqrt(pidn0s * n0sfac / (qs[k] * den[k])));
                rslope[k][1] = wlit<Real>(1.0) / lamdas;
                rslopeb[k][1] = libm::pow(rslope[k][1], bvts);
                rslope2[k][1] = rslope[k][1] * rslope[k][1];
                rslope3[k][1] = rslope2[k][1] * rslope[k][1];
            }

            // Graupel
            if (qg[k] <= qcrmin) {
                rslope[k][2] = rslopegmax;
                rslopeb[k][2] = rslopegbmax;
                rslope2[k][2] = rslopeg2max;
                rslope3[k][2] = rslopeg3max;
            } else {
                Real lamdag = std::sqrt(std::sqrt(pidn0g / (qg[k] * den[k])));
                rslope[k][2] = wlit<Real>(1.0) / lamdag;
                rslopeb[k][2] = libm::pow(rslope[k][2], bvtg);
                rslope2[k][2] = rslope[k][2] * rslope[k][2];
                rslope3[k][2] = rslope2[k][2] * rslope[k][2];
            }

            vt[k][0] = pvtr * rslopeb[k][0] * denfac[k];
            vt[k][1] = pvts * rslopeb[k][1] * denfac[k];
            vt[k][2] = pvtg * rslopeb[k][2] * denfac[k];

            if (qr[k] <= wlit<Real>(0.0)) vt[k][0] = wlit<Real>(0.0);
            if (qs[k] <= wlit<Real>(0.0)) vt[k][1] = wlit<Real>(0.0);
            if (qg[k] <= wlit<Real>(0.0)) vt[k][2] = wlit<Real>(0.0);
        }
    }

    CPPMPAS_HD inline void nislfv_rain_plm(int km, const Real* den, const Real* denfac,
                                        const Real* t, const Real* dz,
                                        Real* ww, Real* rql, Real& precip,
                                        Real dt, int id, int iter,
                                        Real* zi, Real* za, Real* dza, Real* qa,
                                        Real* wi_qmi, Real* wd_qpi, Real* wa, Real* was) const {
        precip = wlit<Real>(0.0);
        Real allold = wlit<Real>(0.0);
        for (int k = 0; k < km; ++k) {
            allold += rql[k];
        }
        if (allold <= wlit<Real>(0.0)) {
            return;
        }

        Real* wi = wi_qmi;
        Real* wd = wd_qpi;

        zi[0] = wlit<Real>(0.0);
        for (int k = 0; k < km; ++k) {
            zi[k + 1] = zi[k] + dz[k];
        }

        for (int k = 0; k < km; ++k) {
            wd[k] = ww[k];
        }

        for (int n = 1; ; ++n) {
            // 2nd order wi
            wi[0] = ww[0];
            wi[km] = ww[km - 1];
            for (int k = 1; k < km; ++k) {
                wi[k] = (ww[k] * dz[k - 1] + ww[k - 1] * dz[k]) / (dz[k - 1] + dz[k]);
            }
            // 3rd order wi
            constexpr Real fa1 = wlit<Real>(9.0) / wlit<Real>(16.0);
            constexpr Real fa2 = wlit<Real>(1.0) / wlit<Real>(16.0);
            wi[0] = ww[0];
            wi[1] = wlit<Real>(0.5) * (ww[1] + ww[0]);
            for (int k = 2; k < km - 1; ++k) {
                wi[k] = fa1 * (ww[k] + ww[k - 1]) - fa2 * (ww[k + 1] + ww[k - 2]);
            }
            if (km >= 2) {
                wi[km - 1] = wlit<Real>(0.5) * (ww[km - 1] + ww[km - 2]);
            }
            wi[km] = ww[km - 1];

            // terminate of top of raingroup
            for (int k = 1; k < km; ++k) {
                if (ww[k] == wlit<Real>(0.0)) {
                    wi[k] = ww[k - 1];
                }
            }

            // diffusivity of wi
            constexpr Real con1 = wlit<Real>(0.05);
            for (int k = km - 1; k >= 0; --k) {
                Real decfl = (wi[k + 1] - wi[k]) * dt / dz[k];
                if (decfl > con1) {
                    wi[k] = wi[k + 1] - con1 * dz[k] / dt;
                }
            }

            // arrival point
            for (int k = 0; k <= km; ++k) {
                za[k] = zi[k] - wi[k] * dt;
            }

            for (int k = 0; k < km; ++k) {
                dza[k] = za[k + 1] - za[k];
            }
            dza[km] = zi[km] - za[km];

            for (int k = 0; k < km; ++k) {
                qa[k] = rql[k] * dz[k] / dza[k];
            }
            qa[km] = wlit<Real>(0.0);

            if (n <= iter) {
                for (int k = 0; k < km; ++k) {
                    Real qr_k = qa[k] / den[k];
                    Real rslope, rslopeb;
                    if (qr_k <= qcrmin) {
                        rslope = rslopermax;
                        rslopeb = rsloperbmax;
                    } else {
                        Real lamdar = std::sqrt(std::sqrt(pidn0r / (qr_k * den[k])));
                        rslope = wlit<Real>(1.0) / lamdar;
                        rslopeb = libm::pow(rslope, bvtr);
                    }
                    wa[k] = pvtr * rslopeb * denfac[k];
                    if (qr_k <= wlit<Real>(0.0)) wa[k] = wlit<Real>(0.0);
                }
                if (n >= 2) {
                    for (int k = 0; k < km; ++k) {
                        wa[k] = wlit<Real>(0.5) * (wa[k] + was[k]);
                    }
                }
                for (int k = 0; k < km; ++k) {
                    ww[k] = wlit<Real>(0.5) * (wd[k] + wa[k]);
                    was[k] = wa[k];
                }
            } else {
                break;
            }
        }

        Real* qmi = wi_qmi;
        Real* qpi = wd_qpi;

        // estimate values at arrival cell interface with monotone
        for (int k = 1; k < km; ++k) {
            Real dip = (qa[k + 1] - qa[k]) / (dza[k + 1] + dza[k]);
            Real dim = (qa[k] - qa[k - 1]) / (dza[k - 1] + dza[k]);
            if (dip * dim <= wlit<Real>(0.0)) {
                qmi[k] = qa[k];
                qpi[k] = qa[k];
            } else {
                qpi[k] = qa[k] + wlit<Real>(0.5) * (dip + dim) * dza[k];
                qmi[k] = wlit<Real>(2.0) * qa[k] - qpi[k];
                if (qpi[k] < wlit<Real>(0.0) || qmi[k] < wlit<Real>(0.0)) {
                    qpi[k] = qa[k];
                    qmi[k] = qa[k];
                }
            }
        }
        qpi[0] = qa[0];
        qmi[0] = qa[0];
        qmi[km] = qa[km];
        qpi[km] = qa[km];

        // interpolation to regular point
        for (int k = 0; k < km; ++k) {
            rql[k] = wlit<Real>(0.0);
        }
        int kb = 0;
        int kt = 0;
        for (int k = 0; k < km; ++k) {
            kb = (kb - 1 > 0) ? kb - 1 : 0;
            kt = (kt - 1 > 0) ? kt - 1 : 0;
            if (zi[k] >= za[km]) {
                break;
            }
            for (int kk = kb; kk < km; ++kk) {
                if (zi[k] <= za[kk + 1]) {
                    kb = kk;
                    break;
                }
            }
            for (int kk = kt; kk < km; ++kk) {
                if (zi[k + 1] <= za[kk]) {
                    kt = kk;
                    break;
                }
            }
            kt = kt - 1;

            if (kt == kb) {
                Real tl = (zi[k] - za[kb]) / dza[kb];
                Real th = (zi[k + 1] - za[kb]) / dza[kb];
                Real tl2 = tl * tl;
                Real th2 = th * th;
                Real qqd = wlit<Real>(0.5) * (qpi[kb] - qmi[kb]);
                Real qqh = qqd * th2 + qmi[kb] * th;
                Real qql = qqd * tl2 + qmi[kb] * tl;
                rql[k] = (qqh - qql) / (th - tl);
            } else if (kt > kb) {
                Real tl = (zi[k] - za[kb]) / dza[kb];
                Real tl2 = tl * tl;
                Real qqd = wlit<Real>(0.5) * (qpi[kb] - qmi[kb]);
                Real qql = qqd * tl2 + qmi[kb] * tl;
                Real dql = qa[kb] - qql;
                Real zsum = (wlit<Real>(1.0) - tl) * dza[kb];
                Real qsum = dql * dza[kb];
                if (kt - kb > 1) {
                    for (int m = kb + 1; m <= kt - 1; ++m) {
                        zsum += dza[m];
                        qsum += qa[m] * dza[m];
                    }
                }
                Real th = (zi[k + 1] - za[kt]) / dza[kt];
                Real th2 = th * th;
                qqd = wlit<Real>(0.5) * (qpi[kt] - qmi[kt]);
                Real dqh = qqd * th2 + qmi[kt] * th;
                zsum += th * dza[kt];
                qsum += dqh * dza[kt];
                rql[k] = qsum / zsum;
            }
        }

        // rain out
        for (int k = 0; k < km; ++k) {
            if (za[k] < wlit<Real>(0.0) && za[k + 1] < wlit<Real>(0.0)) {
                precip += qa[k] * dza[k];
            } else if (za[k] < wlit<Real>(0.0) && za[k + 1] >= wlit<Real>(0.0)) {
                precip += qa[k] * (wlit<Real>(0.0) - za[k]);
                break;
            } else {
                break;
            }
        }
    }

    CPPMPAS_HD inline void nislfv_rain_plm6(int km, const Real* den, const Real* denfac,
                                         const Real* t, const Real* dz,
                                         Real* ww, Real* rql, Real* rql2,
                                         Real& precip1, Real& precip2,
                                         Real dt, int id, int iter,
                                         Real* zi, Real* za, Real* dza, Real* qa, Real* qa2,
                                         Real* wi_qmi, Real* wd_qpi, Real* wa, Real* was) const {
        precip1 = wlit<Real>(0.0);
        precip2 = wlit<Real>(0.0);
        Real allold = wlit<Real>(0.0);
        for (int k = 0; k < km; ++k) {
            allold += rql[k] + rql2[k];
        }
        if (allold <= wlit<Real>(0.0)) {
            return;
        }

        Real* wi = wi_qmi;
        Real* wd = wd_qpi;

        zi[0] = wlit<Real>(0.0);
        for (int k = 0; k < km; ++k) {
            zi[k + 1] = zi[k] + dz[k];
        }

        for (int k = 0; k < km; ++k) {
            wd[k] = ww[k];
        }

        for (int n = 1; ; ++n) {
            // 2nd order wi
            wi[0] = ww[0];
            wi[km] = ww[km - 1];
            for (int k = 1; k < km; ++k) {
                wi[k] = (ww[k] * dz[k - 1] + ww[k - 1] * dz[k]) / (dz[k - 1] + dz[k]);
            }
            // 3rd order wi
            constexpr Real fa1 = wlit<Real>(9.0) / wlit<Real>(16.0);
            constexpr Real fa2 = wlit<Real>(1.0) / wlit<Real>(16.0);
            wi[0] = ww[0];
            wi[1] = wlit<Real>(0.5) * (ww[1] + ww[0]);
            for (int k = 2; k < km - 1; ++k) {
                wi[k] = fa1 * (ww[k] + ww[k - 1]) - fa2 * (ww[k + 1] + ww[k - 2]);
            }
            if (km >= 2) {
                wi[km - 1] = wlit<Real>(0.5) * (ww[km - 1] + ww[km - 2]);
            }
            wi[km] = ww[km - 1];

            // terminate of top of raingroup
            for (int k = 1; k < km; ++k) {
                if (ww[k] == wlit<Real>(0.0)) {
                    wi[k] = ww[k - 1];
                }
            }

            // diffusivity of wi
            constexpr Real con1 = wlit<Real>(0.05);
            for (int k = km - 1; k >= 0; --k) {
                Real decfl = (wi[k + 1] - wi[k]) * dt / dz[k];
                if (decfl > con1) {
                    wi[k] = wi[k + 1] - con1 * dz[k] / dt;
                }
            }

            // arrival point
            for (int k = 0; k <= km; ++k) {
                za[k] = zi[k] - wi[k] * dt;
            }

            for (int k = 0; k < km; ++k) {
                dza[k] = za[k + 1] - za[k];
            }
            dza[km] = zi[km] - za[km];

            for (int k = 0; k < km; ++k) {
                qa[k] = rql[k] * dz[k] / dza[k];
                qa2[k] = rql2[k] * dz[k] / dza[k];
            }
            qa[km] = wlit<Real>(0.0);
            qa2[km] = wlit<Real>(0.0);

            if (n <= iter) {
                for (int k = 0; k < km; ++k) {
                    Real qr_k = qa[k] / den[k];
                    Real qr2_k = qa2[k] / den[k];
                    Real supcol = t0c - t[k];
                    Real n0sfac = std::max(std::min(libm::exp(alpha * supcol), n0smax / n0s), static_cast<Real>(1.0));
                    Real rslope_s, rslopeb_s;
                    if (qr_k <= qcrmin) {
                        rslope_s = rslopesmax;
                        rslopeb_s = rslopesbmax;
                    } else {
                        Real lamdas = std::sqrt(std::sqrt(pidn0s * n0sfac / (qr_k * den[k])));
                        rslope_s = wlit<Real>(1.0) / lamdas;
                        rslopeb_s = libm::pow(rslope_s, bvts);
                    }
                    Real wa_s = pvts * rslopeb_s * denfac[k];
                    if (qr_k <= wlit<Real>(0.0)) wa_s = wlit<Real>(0.0);

                    Real rslope_g, rslopeb_g;
                    if (qr2_k <= qcrmin) {
                        rslope_g = rslopegmax;
                        rslopeb_g = rslopegbmax;
                    } else {
                        Real lamdag = std::sqrt(std::sqrt(pidn0g / (qr2_k * den[k])));
                        rslope_g = wlit<Real>(1.0) / lamdag;
                        rslopeb_g = libm::pow(rslope_g, bvtg);
                    }
                    Real wa_g = pvtg * rslopeb_g * denfac[k];
                    if (qr2_k <= wlit<Real>(0.0)) wa_g = wlit<Real>(0.0);

                    Real tmp_k = std::max(qr_k + qr2_k, static_cast<Real>(1.0e-15));
                    if (tmp_k > static_cast<Real>(1.0e-15)) {
                        wa[k] = (wa_s * qr_k + wa_g * qr2_k) / tmp_k;
                    } else {
                        wa[k] = wlit<Real>(0.0);
                    }
                }
                if (n >= 2) {
                    for (int k = 0; k < km; ++k) {
                        wa[k] = wlit<Real>(0.5) * (wa[k] + was[k]);
                    }
                }
                for (int k = 0; k < km; ++k) {
                    ww[k] = wlit<Real>(0.5) * (wd[k] + wa[k]);
                    was[k] = wa[k];
                }
            } else {
                break;
            }
        }

        Real* qmi = wi_qmi;
        Real* qpi = wd_qpi;

        for (int ist = 1; ist <= 2; ++ist) {
            const Real* cur_qa = (ist == 1) ? qa : qa2;
            Real* out_rql = (ist == 1) ? rql : rql2;
            Real cur_precip = wlit<Real>(0.0);

            for (int k = 1; k < km; ++k) {
                Real dip = (cur_qa[k + 1] - cur_qa[k]) / (dza[k + 1] + dza[k]);
                Real dim = (cur_qa[k] - cur_qa[k - 1]) / (dza[k - 1] + dza[k]);
                if (dip * dim <= wlit<Real>(0.0)) {
                    qmi[k] = cur_qa[k];
                    qpi[k] = cur_qa[k];
                } else {
                    qpi[k] = cur_qa[k] + wlit<Real>(0.5) * (dip + dim) * dza[k];
                    qmi[k] = wlit<Real>(2.0) * cur_qa[k] - qpi[k];
                    if (qpi[k] < wlit<Real>(0.0) || qmi[k] < wlit<Real>(0.0)) {
                        qpi[k] = cur_qa[k];
                        qmi[k] = cur_qa[k];
                    }
                }
            }
            qpi[0] = cur_qa[0];
            qmi[0] = cur_qa[0];
            qmi[km] = cur_qa[km];
            qpi[km] = cur_qa[km];

            for (int k = 0; k < km; ++k) {
                out_rql[k] = wlit<Real>(0.0);
            }
            int kb = 0;
            int kt = 0;
            for (int k = 0; k < km; ++k) {
                kb = (kb - 1 > 0) ? kb - 1 : 0;
                kt = (kt - 1 > 0) ? kt - 1 : 0;
                if (zi[k] >= za[km]) {
                    break;
                }
                for (int kk = kb; kk < km; ++kk) {
                    if (zi[k] <= za[kk + 1]) {
                        kb = kk;
                        break;
                    }
                }
                for (int kk = kt; kk < km; ++kk) {
                    if (zi[k + 1] <= za[kk]) {
                        kt = kk;
                        break;
                    }
                }
                kt = kt - 1;

                if (kt == kb) {
                    Real tl = (zi[k] - za[kb]) / dza[kb];
                    Real th = (zi[k + 1] - za[kb]) / dza[kb];
                    Real tl2 = tl * tl;
                    Real th2 = th * th;
                    Real qqd = wlit<Real>(0.5) * (qpi[kb] - qmi[kb]);
                    Real qqh = qqd * th2 + qmi[kb] * th;
                    Real qql = qqd * tl2 + qmi[kb] * tl;
                    out_rql[k] = (qqh - qql) / (th - tl);
                } else if (kt > kb) {
                    Real tl = (zi[k] - za[kb]) / dza[kb];
                    Real tl2 = tl * tl;
                    Real qqd = wlit<Real>(0.5) * (qpi[kb] - qmi[kb]);
                    Real qql = qqd * tl2 + qmi[kb] * tl;
                    Real dql = cur_qa[kb] - qql;
                    Real zsum = (wlit<Real>(1.0) - tl) * dza[kb];
                    Real qsum = dql * dza[kb];
                    if (kt - kb > 1) {
                        for (int m = kb + 1; m <= kt - 1; ++m) {
                            zsum += dza[m];
                            qsum += cur_qa[m] * dza[m];
                        }
                    }
                    Real th = (zi[k + 1] - za[kt]) / dza[kt];
                    Real th2 = th * th;
                    qqd = wlit<Real>(0.5) * (qpi[kt] - qmi[kt]);
                    Real dqh = qqd * th2 + qmi[kt] * th;
                    zsum += th * dza[kt];
                    qsum += dqh * dza[kt];
                    out_rql[k] = qsum / zsum;
                }
            }

            // rain out
            for (int k = 0; k < km; ++k) {
                if (za[k] < wlit<Real>(0.0) && za[k + 1] < wlit<Real>(0.0)) {
                    cur_precip += cur_qa[k] * dza[k];
                } else if (za[k] < wlit<Real>(0.0) && za[k + 1] >= wlit<Real>(0.0)) {
                    cur_precip += cur_qa[k] * (wlit<Real>(0.0) - za[k]);
                    break;
                } else {
                    break;
                }
            }

            if (ist == 1) {
                precip1 = cur_precip;
            } else {
                precip2 = cur_precip;
            }
        }
    }

    template <int MaxLevels = kMaxLevels>
    CPPMPAS_HD inline void operator()(int nz, Real delt,
                                    Real* th, Real* q, Real* qc, Real* qr,
                                    Real* qi, Real* qs, Real* qg,
                                    const Real* den, const Real* pii, const Real* p,
                                    const Real* delz,
                                    Real& rain, Real& rainncv,
                                    Real& snow, Real& snowncv,
                                    Real& graupel, Real& graupelncv,
                                    Real& sr,
                                    int lane_i = 0, int lane_n = 1) const {  // MPAS parameter (header comment)
        for (int k = 0; k < nz; ++k) {
            qc[k] = std::max(qc[k], static_cast<Real>(0.0));
            qr[k] = std::max(qr[k], static_cast<Real>(0.0));
            qi[k] = std::max(qi[k], static_cast<Real>(0.0));
            qs[k] = std::max(qs[k], static_cast<Real>(0.0));
            qg[k] = std::max(qg[k], static_cast<Real>(0.0));
        }

        Real t[MaxLevels];
        Real cpm[MaxLevels];
        Real xl[MaxLevels];
        for (int k = 0; k < nz; ++k) {
            t[k] = th[k] * pii[k];
            cpm[k] = cpmcal(q[k]);
            xl[k] = xlcal(t[k]);
        }

        rainncv = wlit<Real>(0.0);
        snowncv = wlit<Real>(0.0);
        graupelncv = wlit<Real>(0.0);
        sr = wlit<Real>(0.0);
        Real tstepsnow = wlit<Real>(0.0);
        Real tstepgraup = wlit<Real>(0.0);

        int loops = std::max(1, static_cast<int>(std::round(delt / dtcldcr)));
        Real dtcld = delt / static_cast<Real>(loops);
        if (delt <= dtcldcr) dtcld = delt;

        for (int loop = 0; loop < loops; ++loop) {
            Real denfac[MaxLevels];
            for (int k = 0; k < nz; ++k) {
                denfac[k] = std::sqrt((static_cast<Real>(1.0) / den[k]) * den0);
            }

            Real qsat[MaxLevels][2];
            Real rh[MaxLevels][2];
            Real hsub = xls;
            Real hvap = xlv0;
            Real cvap = cpv;
            Real ttp = t0c + wlit<Real>(0.01);
            Real dldt = cvap - cliq;
            Real xa = -dldt / rv;
            Real xb = xa + hvap / (rv * ttp);
            Real dldti = cvap - cice;
            Real xai = -dldti / rv;
            Real xbi = xai + hsub / (rv * ttp);

            for (int k = 0; k < nz; ++k) {
                Real tr = ttp / t[k];
                Real qs1 = psat * libm::exp(libm::log(tr) * xa) * libm::exp(xb * (wlit<Real>(1.0) - tr));
                qs1 = std::min(qs1, static_cast<Real>(0.99) * p[k]);
                qs1 = ep2 * qs1 / (p[k] - qs1);
                qsat[k][0] = std::max(qs1, qmin);
                rh[k][0] = std::max(q[k] / qsat[k][0], qmin);

                Real qs2;
                if (t[k] < ttp) {
                    qs2 = psat * libm::exp(libm::log(tr) * xai) * libm::exp(xbi * (wlit<Real>(1.0) - tr));
                } else {
                    qs2 = psat * libm::exp(libm::log(tr) * xa) * libm::exp(xb * (wlit<Real>(1.0) - tr));
                }
                qs2 = std::min(qs2, static_cast<Real>(0.99) * p[k]);
                qs2 = ep2 * qs2 / (p[k] - qs2);
                qsat[k][1] = std::max(qs2, qmin);
                rh[k][1] = std::max(q[k] / qsat[k][1], qmin);
            }

            Real xni[MaxLevels];
            for (int k = 0; k < nz; ++k) {
                Real temp = den[k] * std::max(qi[k], qmin);
                temp = std::sqrt(std::sqrt(temp * temp * temp));
                xni[k] = std::min(std::max(static_cast<Real>(5.38e7) * temp, static_cast<Real>(1.0e3)), static_cast<Real>(1.0e6));
            }

            Real rslope[MaxLevels][3];
            Real rslopeb[MaxLevels][3];
            Real rslope2[MaxLevels][3];
            Real rslope3[MaxLevels][3];
            Real vt[MaxLevels][3];
            slope_wsm6(nz, qr, qs, qg, den, denfac, t, rslope, rslopeb, rslope2, rslope3, vt);

            Real workr[MaxLevels];
            Real worka[MaxLevels];
            Real denqrs1[MaxLevels];
            Real denqrs2[MaxLevels];
            Real denqrs3[MaxLevels];
            for (int k = nz - 1; k >= 0; --k) {
                workr[k] = vt[k][0];
                Real qsum_val = std::max(qs[k] + qg[k], static_cast<Real>(1.0e-15));
                if (qsum_val > static_cast<Real>(1.0e-15)) {
                    worka[k] = (vt[k][1] * qs[k] + vt[k][2] * qg[k]) / qsum_val;
                } else {
                    worka[k] = wlit<Real>(0.0);
                }
                denqrs1[k] = den[k] * qr[k];
                denqrs2[k] = den[k] * qs[k];
                denqrs3[k] = den[k] * qg[k];
                if (qr[k] <= wlit<Real>(0.0)) workr[k] = wlit<Real>(0.0);
            }

            Real plm_zi[MaxLevels + 1];
            Real plm_za[MaxLevels + 1];
            Real plm_dza[MaxLevels + 1];
            Real plm_qa[MaxLevels + 1];
            Real plm_qa2[MaxLevels + 1];
            Real plm_wi[MaxLevels + 1];
            Real plm_wd[MaxLevels + 1];
            Real plm_wa[MaxLevels];
            Real plm_was[MaxLevels];

            Real delqrs1 = wlit<Real>(0.0), delqrs2 = wlit<Real>(0.0), delqrs3 = wlit<Real>(0.0);
            nislfv_rain_plm(nz, den, denfac, t, delz, workr, denqrs1, delqrs1, dtcld, 1, 1,
                            plm_zi, plm_za, plm_dza, plm_qa, plm_wi, plm_wd, plm_wa, plm_was);
            nislfv_rain_plm6(nz, den, denfac, t, delz, worka, denqrs2, denqrs3, delqrs2, delqrs3, dtcld, 1, 1,
                             plm_zi, plm_za, plm_dza, plm_qa, plm_qa2, plm_wi, plm_wd, plm_wa, plm_was);

            for (int k = 0; k < nz; ++k) {
                qr[k] = std::max(denqrs1[k] / den[k], static_cast<Real>(0.0));
                qs[k] = std::max(denqrs2[k] / den[k], static_cast<Real>(0.0));
                qg[k] = std::max(denqrs3[k] / den[k], static_cast<Real>(0.0));
            }
            Real fall0_0 = delqrs1 / delz[0] / dtcld;
            Real fall0_1 = delqrs2 / delz[0] / dtcld;
            Real fall0_2 = delqrs3 / delz[0] / dtcld;

            slope_wsm6(nz, qr, qs, qg, den, denfac, t, rslope, rslopeb, rslope2, rslope3, vt);

            for (int k = nz - 1; k >= 0; --k) {
                Real supcol = t0c - t[k];
                Real n0sfac_val = std::max(std::min(libm::exp(alpha * supcol), n0smax / n0s), static_cast<Real>(1.0));
                if (t[k] > t0c) {
                    Real xlf = xlf0;
                    Real work2_val = venfac(p[k], t[k], den[k]);
                    if (qs[k] > wlit<Real>(0.0)) {
                        Real coeres = rslope2[k][1] * std::sqrt(rslope[k][1] * rslopeb[k][1]);
                        Real psmlt_val = xka(t[k], den[k]) / xlf * (t0c - t[k]) * (pi / wlit<Real>(2.0))
                            * n0sfac_val * (precs1 * rslope2[k][1] + precs2 * work2_val * coeres) / den[k];
                        psmlt_val = std::min(std::max(psmlt_val * dtcld, -qs[k]), static_cast<Real>(0.0));
                        qs[k] += psmlt_val;
                        qr[k] -= psmlt_val;
                        t[k] += (xlf / cpm[k]) * psmlt_val;
                    }
                    if (qg[k] > wlit<Real>(0.0)) {
                        Real coeres = rslope2[k][2] * std::sqrt(rslope[k][2] * rslopeb[k][2]);
                        Real pgmlt_val = xka(t[k], den[k]) / xlf * (t0c - t[k])
                            * (precg1 * rslope2[k][2] + precg2 * work2_val * coeres) / den[k];
                        pgmlt_val = std::min(std::max(pgmlt_val * dtcld, -qg[k]), static_cast<Real>(0.0));
                        qg[k] += pgmlt_val;
                        qr[k] -= pgmlt_val;
                        t[k] += (xlf / cpm[k]) * pgmlt_val;
                    }
                }
            }

            Real* work1c = workr;
            Real* denqci = denqrs1;
            for (int k = nz - 1; k >= 0; --k) {
                if (qi[k] <= wlit<Real>(0.0)) {
                    work1c[k] = wlit<Real>(0.0);
                } else {
                    Real xmi = den[k] * qi[k] / xni[k];
                    Real diameter = std::max(std::min(dicon * std::sqrt(xmi), dimax), static_cast<Real>(1.0e-25));
                    work1c[k] = wlit<Real>(1.49e4) * libm::exp(libm::log(diameter) * wlit<Real>(1.31));
                }
                denqci[k] = den[k] * qi[k];
            }
            Real delqi = wlit<Real>(0.0);
            nislfv_rain_plm(nz, den, denfac, t, delz, work1c, denqci, delqi, dtcld, 1, 0,
                            plm_zi, plm_za, plm_dza, plm_qa, plm_wi, plm_wd, plm_wa, plm_was);
            for (int k = 0; k < nz; ++k) {
                qi[k] = std::max(denqci[k] / den[k], static_cast<Real>(0.0));
            }
            Real fallc = delqi / delz[0] / dtcld;

            Real fallsum = fall0_0 + fall0_1 + fall0_2 + fallc;
            Real fallsum_qsi = fall0_1 + fallc;
            Real fallsum_qg = fall0_2;
            if (fallsum > wlit<Real>(0.0)) {
                Real d_rain = fallsum * delz[0] / denr * dtcld * wlit<Real>(1000.0);
                rainncv += d_rain;
                rain += d_rain;
            }
            if (fallsum_qsi > wlit<Real>(0.0)) {
                Real d_snow = fallsum_qsi * delz[0] / denr * dtcld * wlit<Real>(1000.0);
                tstepsnow += d_snow;
                snowncv += d_snow;
                snow += d_snow;
            }
            if (fallsum_qg > wlit<Real>(0.0)) {
                Real d_graup = fallsum_qg * delz[0] / denr * dtcld * wlit<Real>(1000.0);
                tstepgraup += d_graup;
                graupelncv += d_graup;
                graupel += d_graup;
            }
            if (fallsum > wlit<Real>(0.0)) {
                sr = (snowncv + graupelncv) / (rainncv + wlit<Real>(1.0e-12));
            }

            for (int k = 0; k < nz; ++k) {
                Real supcol = t0c - t[k];
                Real xlf = (supcol < wlit<Real>(0.0)) ? xlf0 : (xls - xl[k]);

                // pimlt
                if (supcol < wlit<Real>(0.0) && qi[k] > wlit<Real>(0.0)) {
                    qc[k] += qi[k];
                    t[k] -= (xlf / cpm[k]) * qi[k];
                    qi[k] = wlit<Real>(0.0);
                }
                // pihmf
                if (supcol > wlit<Real>(40.0) && qc[k] > wlit<Real>(0.0)) {
                    qi[k] += qc[k];
                    t[k] += (xlf / cpm[k]) * qc[k];
                    qc[k] = wlit<Real>(0.0);
                }
                // pihtf
                if (supcol > wlit<Real>(0.0) && qc[k] > qmin) {
                    Real supcolt = std::min(supcol, static_cast<Real>(50.0));
                    Real pfrzdtc = std::min(pfrz1 * (libm::exp(pfrz2 * supcolt) - static_cast<Real>(1.0))
                        * den[k] / denr / xncr * qc[k] * qc[k] * dtcld, qc[k]);
                    qi[k] += pfrzdtc;
                    t[k] += (xlf / cpm[k]) * pfrzdtc;
                    qc[k] -= pfrzdtc;
                }
                // pgfrz
                if (supcol > wlit<Real>(0.0) && qr[k] > wlit<Real>(0.0)) {
                    Real temp = rslope3[k][0];
                    temp = temp * temp * rslope[k][0];
                    Real supcolt = std::min(supcol, static_cast<Real>(50.0));
                    Real pfrzdtr = std::min(static_cast<Real>(20.0) * (pi * pi) * pfrz1 * n0r * denr / den[k]
                        * (libm::exp(pfrz2 * supcolt) - static_cast<Real>(1.0)) * temp * dtcld, qr[k]);
                    qg[k] += pfrzdtr;
                    t[k] += (xlf / cpm[k]) * pfrzdtr;
                    qr[k] -= pfrzdtr;
                }
            }

            slope_wsm6(nz, qr, qs, qg, den, denfac, t, rslope, rslopeb, rslope2, rslope3, vt);

            // Warm and cold rain microphysics
            for (int k = 0; k < nz; ++k) {
                // MPAS parameter: the lanes of mp_wsm6.F90:852-854 (diffac / venfac when lane_n = 1)
                const Real work1_0 = diffac_v(xl[k], p[k], t[k], den[k], qsat[k][0], lane_i, lane_n);
                const Real work1_1 = diffac_v(xls, p[k], t[k], den[k], qsat[k][1], lane_i, lane_n);
                const Real work2_k = venfac_v(p[k], t[k], den[k], lane_i, lane_n);
                Real supsat = std::max(q[k], qmin) - qsat[k][0];
                Real satdt = supsat / dtcld;

                Real praut = wlit<Real>(0.0);
                Real pracw = wlit<Real>(0.0);
                Real prevp = wlit<Real>(0.0);

                // praut
                if (qc[k] > qc0) {
                    praut = qck1 * libm::pow(qc[k], static_cast<Real>(7.0 / wlit<Real>(3.0)));
                    praut = std::min(praut, qc[k] / dtcld);
                }
                // pracw
                if (qr[k] > qcrmin && qc[k] > qmin) {
                    pracw = std::min(pacrr * rslope3[k][0] * rslopeb[k][0] * qc[k] * denfac[k], qc[k] / dtcld);
                }
                // prevp
                if (qr[k] > wlit<Real>(0.0)) {
                    Real coeres = rslope2[k][0] * std::sqrt(rslope[k][0] * rslopeb[k][0]);
                    prevp = (rh[k][0] - wlit<Real>(1.0)) * (precr1 * rslope2[k][0] + precr2 * work2_k * coeres) / work1_0;
                    if (prevp < wlit<Real>(0.0)) {
                        prevp = std::max(prevp, -qr[k] / dtcld);
                        prevp = std::max(prevp, satdt / static_cast<Real>(2.0));
                    } else {
                        prevp = std::min(prevp, satdt / static_cast<Real>(2.0));
                    }
                }

                Real supcol = t0c - t[k];
                Real n0sfac_val = std::max(std::min(libm::exp(alpha * supcol), n0smax / n0s), static_cast<Real>(1.0));
                Real supsat_ice = std::max(q[k], qmin) - qsat[k][1];
                Real satdt_ice = supsat_ice / dtcld;
                int ifsat = 0;

                Real temp = den[k] * std::max(qi[k], qmin);
                temp = std::sqrt(std::sqrt(temp * temp * temp));
                xni[k] = std::min(std::max(static_cast<Real>(5.38e7) * temp, static_cast<Real>(1.0e3)), static_cast<Real>(1.0e6));
                Real eacrs = libm::exp(wlit<Real>(0.07) * (-supcol));

                Real xmi = den[k] * qi[k] / xni[k];
                Real diameter = std::min(dicon * std::sqrt(xmi), dimax);
                Real vt2i = wlit<Real>(1.49e4) * libm::pow(diameter, static_cast<Real>(1.31));
                Real vt2r = pvtr * rslopeb[k][0] * denfac[k];
                Real vt2s = pvts * rslopeb[k][1] * denfac[k];
                Real vt2g = pvtg * rslopeb[k][2] * denfac[k];

                Real qsum_val = std::max(qs[k] + qg[k], static_cast<Real>(1.0e-15));
                Real vt2ave = (qsum_val > static_cast<Real>(1.0e-15))
                    ? (vt2s * qs[k] + vt2g * qg[k]) / qsum_val : static_cast<Real>(0.0);

                Real praci = wlit<Real>(0.0);
                Real piacr = wlit<Real>(0.0);
                Real psaci = wlit<Real>(0.0);
                Real pgaci = wlit<Real>(0.0);

                if (supcol > wlit<Real>(0.0) && qi[k] > qmin) {
                    if (qr[k] > qcrmin) {
                        Real acrfac = wlit<Real>(2.0) * rslope3[k][0] + wlit<Real>(2.0) * diameter * rslope2[k][0]
                            + diameter * diameter * rslope[k][0];
                        praci = pi * qi[k] * n0r * std::abs(vt2r - vt2i) * acrfac / wlit<Real>(4.0);
                        Real eff = std::min(std::max(static_cast<Real>(0.0), qr[k] / qi[k]), static_cast<Real>(1.0));
                        praci = praci * (eff * eff);
                        praci = std::min(praci, qi[k] / dtcld);

                        piacr = (pi * pi) * avtr * n0r * denr * xni[k] * denfac[k]
                            * g6pbr * rslope3[k][0] * rslope3[k][0]
                            * rslopeb[k][0] / wlit<Real>(24.0) / den[k];
                        Real eff2 = std::min(std::max(static_cast<Real>(0.0), qi[k] / qr[k]), static_cast<Real>(1.0));
                        piacr = piacr * (eff2 * eff2);
                        piacr = std::min(piacr, qr[k] / dtcld);
                    }
                    if (qs[k] > qcrmin) {
                        Real acrfac = wlit<Real>(2.0) * rslope3[k][1] + wlit<Real>(2.0) * diameter * rslope2[k][1]
                            + diameter * diameter * rslope[k][1];
                        psaci = pi * qi[k] * eacrs * n0s * n0sfac_val * std::abs(vt2ave - vt2i) * acrfac / wlit<Real>(4.0);
                        psaci = std::min(psaci, qi[k] / dtcld);
                    }
                    if (qg[k] > qcrmin) {
                        Real egi = libm::exp(wlit<Real>(0.07) * (-supcol));
                        Real acrfac = wlit<Real>(2.0) * rslope3[k][2] + wlit<Real>(2.0) * diameter * rslope2[k][2]
                            + diameter * diameter * rslope[k][2];
                        pgaci = pi * egi * qi[k] * n0g * std::abs(vt2ave - vt2i) * acrfac / wlit<Real>(4.0);
                        pgaci = std::min(pgaci, qi[k] / dtcld);
                    }
                }

                Real psacw = wlit<Real>(0.0);
                if (qs[k] > qcrmin && qc[k] > qmin) {
                    Real eff = std::min(std::max(static_cast<Real>(0.0), qs[k] / qc[k]), static_cast<Real>(1.0));
                    psacw = std::min(pacrc * n0sfac_val * rslope3[k][1] * rslopeb[k][1]
                        * (eff * eff) * qc[k] * denfac[k], qc[k] / dtcld);
                }

                Real pgacw = wlit<Real>(0.0);
                if (qg[k] > qcrmin && qc[k] > qmin) {
                    Real eff = std::min(std::max(static_cast<Real>(0.0), qg[k] / qc[k]), static_cast<Real>(1.0));
                    pgacw = std::min(pacrg * rslope3[k][2] * rslopeb[k][2]
                        * (eff * eff) * qc[k] * denfac[k], qc[k] / dtcld);
                }

                Real paacw = wlit<Real>(0.0);
                if (qsum_val > static_cast<Real>(1.0e-15)) {
                    paacw = (qs[k] * psacw + qg[k] * pgacw) / qsum_val;
                }

                Real pracs = wlit<Real>(0.0);
                Real psacr = wlit<Real>(0.0);
                if (qs[k] > qcrmin && qr[k] > qcrmin) {
                    if (supcol > wlit<Real>(0.0)) {
                        Real acrfac = wlit<Real>(5.0) * rslope3[k][1] * rslope3[k][1] * rslope[k][0]
                            + wlit<Real>(2.0) * rslope3[k][1] * rslope2[k][1] * rslope2[k][0]
                            + wlit<Real>(0.5) * rslope2[k][1] * rslope2[k][1] * rslope3[k][0];
                        pracs = (pi * pi) * n0r * n0s * n0sfac_val * std::abs(vt2r - vt2ave)
                            * (dens / den[k]) * acrfac;
                        Real eff = std::min(std::max(static_cast<Real>(0.0), qr[k] / qs[k]), static_cast<Real>(1.0));
                        pracs = pracs * (eff * eff);
                        pracs = std::min(pracs, qs[k] / dtcld);
                    }
                    Real acrfac = wlit<Real>(5.0) * rslope3[k][0] * rslope3[k][0] * rslope[k][1]
                        + wlit<Real>(2.0) * rslope3[k][0] * rslope2[k][0] * rslope2[k][1]
                        + wlit<Real>(0.5) * rslope2[k][0] * rslope2[k][0] * rslope3[k][1];
                    psacr = (pi * pi) * n0r * n0s * n0sfac_val * std::abs(vt2ave - vt2r)
                        * (denr / den[k]) * acrfac;
                    Real eff = std::min(std::max(static_cast<Real>(0.0), qs[k] / qr[k]), static_cast<Real>(1.0));
                    psacr = psacr * (eff * eff);
                    psacr = std::min(psacr, qr[k] / dtcld);
                }

                Real pgacr = wlit<Real>(0.0);
                if (qg[k] > qcrmin && qr[k] > qcrmin) {
                    Real acrfac = wlit<Real>(5.0) * rslope3[k][0] * rslope3[k][0] * rslope[k][2]
                        + wlit<Real>(2.0) * rslope3[k][0] * rslope2[k][0] * rslope2[k][2]
                        + wlit<Real>(0.5) * rslope2[k][0] * rslope2[k][0] * rslope3[k][2];
                    pgacr = (pi * pi) * n0r * n0g * std::abs(vt2ave - vt2r) * (denr / den[k]) * acrfac;
                    Real eff = std::min(std::max(static_cast<Real>(0.0), qg[k] / qr[k]), static_cast<Real>(1.0));
                    pgacr = pgacr * (eff * eff);
                    pgacr = std::min(pgacr, qr[k] / dtcld);
                }

                Real pgacs = wlit<Real>(0.0);
                Real pseml = wlit<Real>(0.0);
                Real pgeml = wlit<Real>(0.0);
                if (supcol <= wlit<Real>(0.0)) {
                    Real xlf = xlf0;
                    if (qs[k] > wlit<Real>(0.0)) {
                        pseml = std::min(std::max(cliq * supcol * (paacw + psacr) / xlf, -qs[k] / dtcld), static_cast<Real>(0.0));
                    }
                    if (qg[k] > wlit<Real>(0.0)) {
                        pgeml = std::min(std::max(cliq * supcol * (paacw + pgacr) / xlf, -qg[k] / dtcld), static_cast<Real>(0.0));
                    }
                }

                Real pidep = wlit<Real>(0.0);
                Real psdep = wlit<Real>(0.0);
                Real pgdep = wlit<Real>(0.0);
                Real pigen = wlit<Real>(0.0);
                Real psaut = wlit<Real>(0.0);
                Real pgaut = wlit<Real>(0.0);

                if (supcol > wlit<Real>(0.0)) {
                    if (qi[k] > wlit<Real>(0.0) && ifsat != 1) {
                        pidep = wlit<Real>(4.0) * diameter * xni[k] * (rh[k][1] - wlit<Real>(1.0)) / work1_1;
                        Real supice = satdt_ice - prevp;
                        if (pidep < wlit<Real>(0.0)) {
                            pidep = std::max(std::max(pidep, satdt_ice / static_cast<Real>(2.0)), supice);
                            pidep = std::max(pidep, -qi[k] / dtcld);
                        } else {
                            pidep = std::min(std::min(pidep, satdt_ice / static_cast<Real>(2.0)), supice);
                        }
                        if (std::abs(prevp + pidep) >= std::abs(satdt_ice)) ifsat = 1;
                    }

                    if (qs[k] > wlit<Real>(0.0) && ifsat != 1) {
                        Real coeres = rslope2[k][1] * std::sqrt(rslope[k][1] * rslopeb[k][1]);
                        psdep = (rh[k][1] - wlit<Real>(1.0)) * n0sfac_val * (precs1 * rslope2[k][1]
                            + precs2 * work2_k * coeres) / work1_1;
                        Real supice = satdt_ice - prevp - pidep;
                        if (psdep < wlit<Real>(0.0)) {
                            psdep = std::max(psdep, -qs[k] / dtcld);
                            psdep = std::max(std::max(psdep, satdt_ice / static_cast<Real>(2.0)), supice);
                        } else {
                            psdep = std::min(std::min(psdep, satdt_ice / static_cast<Real>(2.0)), supice);
                        }
                        if (std::abs(prevp + pidep + psdep) >= std::abs(satdt_ice)) ifsat = 1;
                    }

                    if (qg[k] > wlit<Real>(0.0) && ifsat != 1) {
                        Real coeres = rslope2[k][2] * std::sqrt(rslope[k][2] * rslopeb[k][2]);
                        pgdep = (rh[k][1] - wlit<Real>(1.0)) * (precg1 * rslope2[k][2]
                            + precg2 * work2_k * coeres) / work1_1;
                        Real supice = satdt_ice - prevp - pidep - psdep;
                        if (pgdep < wlit<Real>(0.0)) {
                            pgdep = std::max(pgdep, -qg[k] / dtcld);
                            pgdep = std::max(std::max(pgdep, satdt_ice / static_cast<Real>(2.0)), supice);
                        } else {
                            pgdep = std::min(std::min(pgdep, satdt_ice / static_cast<Real>(2.0)), supice);
                        }
                        if (std::abs(prevp + pidep + psdep + pgdep) >= std::abs(satdt_ice)) ifsat = 1;
                    }

                    if (supsat_ice > wlit<Real>(0.0) && ifsat != 1) {
                        Real supice = satdt_ice - prevp - pidep - psdep - pgdep;
                        Real xni0 = wlit<Real>(1.0e3) * libm::exp(wlit<Real>(0.1) * supcol);
                        Real roqi0 = wlit<Real>(4.92e-11) * libm::pow(xni0, static_cast<Real>(1.33));
                        pigen = std::max(static_cast<Real>(0.0), (roqi0 / den[k] - std::max(qi[k], static_cast<Real>(0.0))) / dtcld);
                        pigen = std::min(std::min(pigen, satdt_ice), supice);
                    }

                    if (qi[k] > wlit<Real>(0.0)) {
                        Real qimax = roqimax / den[k];
                        psaut = std::max(static_cast<Real>(0.0), (qi[k] - qimax) / dtcld);
                    }

                    if (qs[k] > wlit<Real>(0.0)) {
                        Real alpha2 = wlit<Real>(1.0e-3) * libm::exp(wlit<Real>(0.09) * (-supcol));
                        pgaut = std::min(std::max(static_cast<Real>(0.0), alpha2 * (qs[k] - qs0)), qs[k] / dtcld);
                    }
                }

                Real psevp = wlit<Real>(0.0);
                Real pgevp = wlit<Real>(0.0);
                if (supcol < wlit<Real>(0.0)) {
                    if (qs[k] > wlit<Real>(0.0) && rh[k][0] < wlit<Real>(1.0)) {
                        Real coeres = rslope2[k][1] * std::sqrt(rslope[k][1] * rslopeb[k][1]);
                        psevp = (rh[k][0] - wlit<Real>(1.0)) * n0sfac_val * (precs1 * rslope2[k][1]
                            + precs2 * work2_k * coeres) / work1_0;
                        psevp = std::min(std::max(psevp, -qs[k] / dtcld), static_cast<Real>(0.0));
                    }
                    if (qg[k] > wlit<Real>(0.0) && rh[k][0] < wlit<Real>(1.0)) {
                        Real coeres = rslope2[k][2] * std::sqrt(rslope[k][2] * rslopeb[k][2]);
                        pgevp = (rh[k][0] - wlit<Real>(1.0)) * (precg1 * rslope2[k][2]
                            + precg2 * work2_k * coeres) / work1_0;
                        pgevp = std::min(std::max(pgevp, -qg[k] / dtcld), static_cast<Real>(0.0));
                    }
                }

                // Mass conservation check and feedback
                Real delta2 = (qr[k] < wlit<Real>(1.0e-4) && qs[k] < wlit<Real>(1.0e-4)) ? static_cast<Real>(1.0) : static_cast<Real>(0.0);
                Real delta3 = (qr[k] < wlit<Real>(1.0e-4)) ? static_cast<Real>(1.0) : static_cast<Real>(0.0);

                if (t[k] <= t0c) {
                    // cloud water
                    Real value = std::max(qmin, qc[k]);
                    Real source = (praut + pracw + paacw + paacw) * dtcld;
                    if (source > value) {
                        Real factor = value / source;
                        praut *= factor;
                        pracw *= factor;
                        paacw *= factor;
                    }

                    // cloud ice
                    value = std::max(qmin, qi[k]);
                    source = (psaut - pigen - pidep + praci + psaci + pgaci) * dtcld;
                    if (source > value) {
                        Real factor = value / source;
                        psaut *= factor;
                        pigen *= factor;
                        pidep *= factor;
                        praci *= factor;
                        psaci *= factor;
                        pgaci *= factor;
                    }

                    // rain
                    value = std::max(qmin, qr[k]);
                    source = (-praut - prevp - pracw + piacr + psacr + pgacr) * dtcld;
                    if (source > value) {
                        Real factor = value / source;
                        praut *= factor;
                        prevp *= factor;
                        pracw *= factor;
                        piacr *= factor;
                        psacr *= factor;
                        pgacr *= factor;
                    }

                    // snow
                    value = std::max(qmin, qs[k]);
                    source = -(psdep + psaut - pgaut + paacw + piacr * delta3
                        + praci * delta3 - pracs * (wlit<Real>(1.0) - delta2)
                        + psacr * delta2 + psaci - pgacs) * dtcld;
                    if (source > value) {
                        Real factor = value / source;
                        psdep *= factor;
                        psaut *= factor;
                        pgaut *= factor;
                        paacw *= factor;
                        piacr *= factor;
                        praci *= factor;
                        psaci *= factor;
                        pracs *= factor;
                        psacr *= factor;
                        pgacs *= factor;
                    }

                    // graupel
                    value = std::max(qmin, qg[k]);
                    source = -(pgdep + pgaut + piacr * (wlit<Real>(1.0) - delta3) + praci * (wlit<Real>(1.0) - delta3)
                        + psacr * (wlit<Real>(1.0) - delta2) + pracs * (wlit<Real>(1.0) - delta2)
                        + pgaci + paacw + pgacr + pgacs) * dtcld;
                    if (source > value) {
                        Real factor = value / source;
                        pgdep *= factor;
                        pgaut *= factor;
                        piacr *= factor;
                        praci *= factor;
                        psacr *= factor;
                        pracs *= factor;
                        paacw *= factor;
                        pgaci *= factor;
                        pgacr *= factor;
                        pgacs *= factor;
                    }

                    Real work2_val = -(prevp + psdep + pgdep + pigen + pidep);
                    q[k] += work2_val * dtcld;
                    qc[k] = std::max(qc[k] - (praut + pracw + paacw + paacw) * dtcld, static_cast<Real>(0.0));
                    qr[k] = std::max(qr[k] + (praut + pracw + prevp - piacr - pgacr - psacr) * dtcld, static_cast<Real>(0.0));
                    qi[k] = std::max(qi[k] - (psaut + praci + psaci + pgaci - pigen - pidep) * dtcld, static_cast<Real>(0.0));
                    qs[k] = std::max(qs[k] + (psdep + psaut + paacw - pgaut + piacr * delta3
                        + praci * delta3 + psaci - pgacs - pracs * (static_cast<Real>(1.0) - delta2) + psacr * delta2) * dtcld, static_cast<Real>(0.0));
                    qg[k] = std::max(qg[k] + (pgdep + pgaut + piacr * (static_cast<Real>(1.0) - delta3)
                        + praci * (static_cast<Real>(1.0) - delta3) + psacr * (static_cast<Real>(1.0) - delta2)
                        + pracs * (static_cast<Real>(1.0) - delta2) + pgaci + paacw + pgacr + pgacs) * dtcld, static_cast<Real>(0.0));

                    Real xlf = xls - xl[k];
                    Real xlwork2 = -xls * (psdep + pgdep + pidep + pigen)
                        - xl[k] * prevp - xlf * (piacr + paacw + paacw + pgacr + psacr);
                    t[k] -= xlwork2 / cpm[k] * dtcld;

                } else {
                    // t[k] > t0c
                    // cloud water
                    Real value = std::max(qmin, qc[k]);
                    Real source = (praut + pracw + paacw + paacw) * dtcld;
                    if (source > value) {
                        Real factor = value / source;
                        praut *= factor;
                        pracw *= factor;
                        paacw *= factor;
                    }

                    // rain
                    value = std::max(qmin, qr[k]);
                    source = (-paacw - praut + pseml + pgeml - pracw - paacw - prevp) * dtcld;
                    if (source > value) {
                        Real factor = value / source;
                        praut *= factor;
                        prevp *= factor;
                        pracw *= factor;
                        paacw *= factor;
                        pseml *= factor;
                        pgeml *= factor;
                    }

                    // snow
                    value = std::max(qcrmin, qs[k]);
                    source = (pgacs - pseml - psevp) * dtcld;
                    if (source > value) {
                        Real factor = value / source;
                        pgacs *= factor;
                        psevp *= factor;
                        pseml *= factor;
                    }

                    // graupel
                    value = std::max(qcrmin, qg[k]);
                    source = -(pgacs + pgevp + pgeml) * dtcld;
                    if (source > value) {
                        Real factor = value / source;
                        pgacs *= factor;
                        pgevp *= factor;
                        pgeml *= factor;
                    }

                    Real work2_val = -(prevp + psevp + pgevp);
                    q[k] += work2_val * dtcld;
                    qc[k] = std::max(qc[k] - (praut + pracw + paacw + paacw) * dtcld, static_cast<Real>(0.0));
                    qr[k] = std::max(qr[k] + (praut + pracw + prevp + paacw + paacw - pseml - pgeml) * dtcld, static_cast<Real>(0.0));
                    qs[k] = std::max(qs[k] + (psevp - pgacs + pseml) * dtcld, static_cast<Real>(0.0));
                    qg[k] = std::max(qg[k] + (pgacs + pgevp + pgeml) * dtcld, static_cast<Real>(0.0));

                    Real xlf = xls - xl[k];
                    Real xlwork2 = -xl[k] * (prevp + psevp + pgevp) - xlf * (pseml + pgeml);
                    t[k] -= xlwork2 / cpm[k] * dtcld;
                }
            }

            // Saturation adjustment (conden)
            for (int k = 0; k < nz; ++k) {
                Real tr = ttp / t[k];
                Real qs1 = psat * libm::exp(libm::log(tr) * xa) * libm::exp(xb * (wlit<Real>(1.0) - tr));
                qs1 = std::min(qs1, static_cast<Real>(0.99) * p[k]);
                qs1 = ep2 * qs1 / (p[k] - qs1);
                qs1 = std::max(qs1, qmin);

                Real work1_cond = conden(t[k], q[k], qs1, xl[k], cpm[k]);
                Real pcond = std::min(std::max(work1_cond / dtcld, static_cast<Real>(0.0)),
                                      std::max(q[k], static_cast<Real>(0.0)) / dtcld);
                if (qc[k] > wlit<Real>(0.0) && work1_cond < wlit<Real>(0.0)) {
                    pcond = std::max(work1_cond, -qc[k]) / dtcld;
                }
                q[k] -= pcond * dtcld;
                qc[k] = std::max(qc[k] + pcond * dtcld, static_cast<Real>(0.0));
                t[k] += pcond * xl[k] / cpm[k] * dtcld;
            }

            // Small value padding
            for (int k = 0; k < nz; ++k) {
                if (qc[k] <= qmin) qc[k] = wlit<Real>(0.0);
                if (qi[k] <= qmin) qi[k] = wlit<Real>(0.0);
            }
        } // end loops

        for (int k = 0; k < nz; ++k) {
            th[k] = t[k] / pii[k];
        }
    }

#if defined(CPPWRF_DOUBLE)
    // ---- WSM6 in column stages (FP64 builds, CPU and CUDA; Wsm6Stages below) ----
    // Each function is the body of one of operator()'s level loops for one level k: the same statements with the
    // column arrays as scalars (x[k] -> x), so every value is the same arithmetic in the same order.

    // the body of slope_wsm6's level loop
    CPPMPAS_HD inline void slope_level(Real qr, Real qs, Real qg, Real den, Real denfac, Real t,
                                      Real rslope[3], Real rslopeb[3], Real rslope2[3], Real rslope3[3],
                                      Real vt[3]) const {
        Real supcol = t0c - t;
        Real n0sfac = std::max(std::min(libm::exp(alpha * supcol), n0smax / n0s), static_cast<Real>(1.0));

        // Rain
        if (qr <= qcrmin) {
            rslope[0] = rslopermax;
            rslopeb[0] = rsloperbmax;
            rslope2[0] = rsloper2max;
            rslope3[0] = rsloper3max;
        } else {
            Real lamdar = std::sqrt(std::sqrt(pidn0r / (qr * den)));
            rslope[0] = wlit<Real>(1.0) / lamdar;
            rslopeb[0] = libm::pow(rslope[0], bvtr);
            rslope2[0] = rslope[0] * rslope[0];
            rslope3[0] = rslope2[0] * rslope[0];
        }

        // Snow
        if (qs <= qcrmin) {
            rslope[1] = rslopesmax;
            rslopeb[1] = rslopesbmax;
            rslope2[1] = rslopes2max;
            rslope3[1] = rslopes3max;
        } else {
            Real lamdas = std::sqrt(std::sqrt(pidn0s * n0sfac / (qs * den)));
            rslope[1] = wlit<Real>(1.0) / lamdas;
            rslopeb[1] = libm::pow(rslope[1], bvts);
            rslope2[1] = rslope[1] * rslope[1];
            rslope3[1] = rslope2[1] * rslope[1];
        }

        // Graupel
        if (qg <= qcrmin) {
            rslope[2] = rslopegmax;
            rslopeb[2] = rslopegbmax;
            rslope2[2] = rslopeg2max;
            rslope3[2] = rslopeg3max;
        } else {
            Real lamdag = std::sqrt(std::sqrt(pidn0g / (qg * den)));
            rslope[2] = wlit<Real>(1.0) / lamdag;
            rslopeb[2] = libm::pow(rslope[2], bvtg);
            rslope2[2] = rslope[2] * rslope[2];
            rslope3[2] = rslope2[2] * rslope[2];
        }

        vt[0] = pvtr * rslopeb[0] * denfac;
        vt[1] = pvts * rslopeb[1] * denfac;
        vt[2] = pvtg * rslopeb[2] * denfac;

        if (qr <= wlit<Real>(0.0)) vt[0] = wlit<Real>(0.0);
        if (qs <= wlit<Real>(0.0)) vt[1] = wlit<Real>(0.0);
        if (qg <= wlit<Real>(0.0)) vt[2] = wlit<Real>(0.0);
    }

    // before the loops: the clamps of qc ... qg, t, cpm, xl
    CPPMPAS_HD inline void level_init(Real th, Real pii, Real q, Real& qc, Real& qr, Real& qi, Real& qs, Real& qg,
                                     Real& t, Real& cpm, Real& xl) const {
        qc = std::max(qc, static_cast<Real>(0.0));
        qr = std::max(qr, static_cast<Real>(0.0));
        qi = std::max(qi, static_cast<Real>(0.0));
        qs = std::max(qs, static_cast<Real>(0.0));
        qg = std::max(qg, static_cast<Real>(0.0));
        t = th * pii;
        cpm = cpmcal(q);
        xl = xlcal(t);
    }

    // start of a loop: qsat, rh; the fall speeds and densities of the rain and snow + graupel sedimentation
    CPPMPAS_HD inline void level_a(Real t, Real q, Real qr, Real qs, Real qg, Real den, Real p, Real qsat[2],
                                  Real rh[2], Real& workr, Real& worka, Real& denqrs1, Real& denqrs2,
                                  Real& denqrs3) const {
        Real denfac;
        denfac = std::sqrt((static_cast<Real>(1.0) / den) * den0);
        Real hsub = xls;
        Real hvap = xlv0;
        Real cvap = cpv;
        Real ttp = t0c + wlit<Real>(0.01);
        Real dldt = cvap - cliq;
        Real xa = -dldt / rv;
        Real xb = xa + hvap / (rv * ttp);
        Real dldti = cvap - cice;
        Real xai = -dldti / rv;
        Real xbi = xai + hsub / (rv * ttp);
        {
            Real tr = ttp / t;
            Real qs1 = psat * libm::exp(libm::log(tr) * xa) * libm::exp(xb * (wlit<Real>(1.0) - tr));
            qs1 = std::min(qs1, static_cast<Real>(0.99) * p);
            qs1 = ep2 * qs1 / (p - qs1);
            qsat[0] = std::max(qs1, qmin);
            rh[0] = std::max(q / qsat[0], qmin);

            Real qs2;
            if (t < ttp) {
                qs2 = psat * libm::exp(libm::log(tr) * xai) * libm::exp(xbi * (wlit<Real>(1.0) - tr));
            } else {
                qs2 = psat * libm::exp(libm::log(tr) * xa) * libm::exp(xb * (wlit<Real>(1.0) - tr));
            }
            qs2 = std::min(qs2, static_cast<Real>(0.99) * p);
            qs2 = ep2 * qs2 / (p - qs2);
            qsat[1] = std::max(qs2, qmin);
            rh[1] = std::max(q / qsat[1], qmin);
        }
        Real rslope[3], rslopeb[3], rslope2[3], rslope3[3], vt[3];
        slope_level(qr, qs, qg, den, denfac, t, rslope, rslopeb, rslope2, rslope3, vt);
        {
            workr = vt[0];
            Real qsum_val = std::max(qs + qg, static_cast<Real>(1.0e-15));
            if (qsum_val > static_cast<Real>(1.0e-15)) {
                worka = (vt[1] * qs + vt[2] * qg) / qsum_val;
            } else {
                worka = wlit<Real>(0.0);
            }
            denqrs1 = den * qr;
            denqrs2 = den * qs;
            denqrs3 = den * qg;
            if (qr <= wlit<Real>(0.0)) workr = wlit<Real>(0.0);
        }
    }

    // after the rain / snow + graupel sedimentation: qr, qs, qg, the melting, the ice fall speed (xni of the
    // loop's start: qi is unchanged since); rs0 / rs30: rslope / rslope3 of rain for the freezing in level_c
    CPPMPAS_HD inline void level_b(Real& t, Real cpm, Real& qr, Real& qs, Real& qg, Real qi, Real den, Real p,
                                  Real denqrs1, Real denqrs2, Real denqrs3, Real dtcld, Real& work1c,
                                  Real& denqci, Real& rs0, Real& rs30) const {
        Real denfac;
        denfac = std::sqrt((static_cast<Real>(1.0) / den) * den0);
        qr = std::max(denqrs1 / den, static_cast<Real>(0.0));
        qs = std::max(denqrs2 / den, static_cast<Real>(0.0));
        qg = std::max(denqrs3 / den, static_cast<Real>(0.0));
        Real rslope[3], rslopeb[3], rslope2[3], rslope3[3], vt[3];
        slope_level(qr, qs, qg, den, denfac, t, rslope, rslopeb, rslope2, rslope3, vt);
        {
            Real supcol = t0c - t;
            Real n0sfac_val = std::max(std::min(libm::exp(alpha * supcol), n0smax / n0s), static_cast<Real>(1.0));
            if (t > t0c) {
                Real xlf = xlf0;
                Real work2_val = venfac(p, t, den);
                if (qs > wlit<Real>(0.0)) {
                    Real coeres = rslope2[1] * std::sqrt(rslope[1] * rslopeb[1]);
                    Real psmlt_val = xka(t, den) / xlf * (t0c - t) * (pi / wlit<Real>(2.0))
                        * n0sfac_val * (precs1 * rslope2[1] + precs2 * work2_val * coeres) / den;
                    psmlt_val = std::min(std::max(psmlt_val * dtcld, -qs), static_cast<Real>(0.0));
                    qs += psmlt_val;
                    qr -= psmlt_val;
                    t += (xlf / cpm) * psmlt_val;
                }
                if (qg > wlit<Real>(0.0)) {
                    Real coeres = rslope2[2] * std::sqrt(rslope[2] * rslopeb[2]);
                    Real pgmlt_val = xka(t, den) / xlf * (t0c - t)
                        * (precg1 * rslope2[2] + precg2 * work2_val * coeres) / den;
                    pgmlt_val = std::min(std::max(pgmlt_val * dtcld, -qg), static_cast<Real>(0.0));
                    qg += pgmlt_val;
                    qr -= pgmlt_val;
                    t += (xlf / cpm) * pgmlt_val;
                }
            }
        }
        Real xni;
        {
            Real temp = den * std::max(qi, qmin);
            temp = std::sqrt(std::sqrt(temp * temp * temp));
            xni = std::min(std::max(static_cast<Real>(5.38e7) * temp, static_cast<Real>(1.0e3)), static_cast<Real>(1.0e6));
        }
        {
            if (qi <= wlit<Real>(0.0)) {
                work1c = wlit<Real>(0.0);
            } else {
                Real xmi = den * qi / xni;
                Real diameter = std::max(std::min(dicon * std::sqrt(xmi), dimax), static_cast<Real>(1.0e-25));
                work1c = wlit<Real>(1.49e4) * libm::exp(libm::log(diameter) * wlit<Real>(1.31));
            }
            denqci = den * qi;
        }
        rs0 = rslope[0];
        rs30 = rslope3[0];
    }

    // after the ice sedimentation: qi, the freezing (rslope / rslope3 of rain from level_b), the warm and cold
    // rain processes, the saturation adjustment, the small-value padding
    CPPMPAS_HD inline void level_c(Real& t, Real& q, Real& qc, Real& qr, Real& qi, Real& qs, Real& qg, Real den,
                                  Real p, Real cpm, Real xl, const Real qsat[2], const Real rh[2], Real denqci,
                                  Real rs0, Real rs30, Real dtcld) const {
        Real denfac;
        denfac = std::sqrt((static_cast<Real>(1.0) / den) * den0);
        qi = std::max(denqci / den, static_cast<Real>(0.0));
        Real rslope[3], rslopeb[3], rslope2[3], rslope3[3], vt[3];
        rslope[0] = rs0;
        rslope3[0] = rs30;
        {
            Real supcol = t0c - t;
            Real xlf = (supcol < wlit<Real>(0.0)) ? xlf0 : (xls - xl);

            // pimlt
            if (supcol < wlit<Real>(0.0) && qi > wlit<Real>(0.0)) {
                qc += qi;
                t -= (xlf / cpm) * qi;
                qi = wlit<Real>(0.0);
            }
            // pihmf
            if (supcol > wlit<Real>(40.0) && qc > wlit<Real>(0.0)) {
                qi += qc;
                t += (xlf / cpm) * qc;
                qc = wlit<Real>(0.0);
            }
            // pihtf
            if (supcol > wlit<Real>(0.0) && qc > qmin) {
                Real supcolt = std::min(supcol, static_cast<Real>(50.0));
                Real pfrzdtc = std::min(pfrz1 * (libm::exp(pfrz2 * supcolt) - static_cast<Real>(1.0))
                    * den / denr / xncr * qc * qc * dtcld, qc);
                qi += pfrzdtc;
                t += (xlf / cpm) * pfrzdtc;
                qc -= pfrzdtc;
            }
            // pgfrz
            if (supcol > wlit<Real>(0.0) && qr > wlit<Real>(0.0)) {
                Real temp = rslope3[0];
                temp = temp * temp * rslope[0];
                Real supcolt = std::min(supcol, static_cast<Real>(50.0));
                Real pfrzdtr = std::min(static_cast<Real>(20.0) * (pi * pi) * pfrz1 * n0r * denr / den
                    * (libm::exp(pfrz2 * supcolt) - static_cast<Real>(1.0)) * temp * dtcld, qr);
                qg += pfrzdtr;
                t += (xlf / cpm) * pfrzdtr;
                qr -= pfrzdtr;
            }
        }
        slope_level(qr, qs, qg, den, denfac, t, rslope, rslopeb, rslope2, rslope3, vt);
        Real xni;
        {
            const Real work1_0 = diffac(xl, p, t, den, qsat[0]);
            const Real work1_1 = diffac(xls, p, t, den, qsat[1]);
            const Real work2_k = venfac(p, t, den);
            Real supsat = std::max(q, qmin) - qsat[0];
            Real satdt = supsat / dtcld;

            Real praut = wlit<Real>(0.0);
            Real pracw = wlit<Real>(0.0);
            Real prevp = wlit<Real>(0.0);

            // praut
            if (qc > qc0) {
                praut = qck1 * libm::pow(qc, static_cast<Real>(7.0 / wlit<Real>(3.0)));
                praut = std::min(praut, qc / dtcld);
            }
            // pracw
            if (qr > qcrmin && qc > qmin) {
                pracw = std::min(pacrr * rslope3[0] * rslopeb[0] * qc * denfac, qc / dtcld);
            }
            // prevp
            if (qr > wlit<Real>(0.0)) {
                Real coeres = rslope2[0] * std::sqrt(rslope[0] * rslopeb[0]);
                prevp = (rh[0] - wlit<Real>(1.0)) * (precr1 * rslope2[0] + precr2 * work2_k * coeres) / work1_0;
                if (prevp < wlit<Real>(0.0)) {
                    prevp = std::max(prevp, -qr / dtcld);
                    prevp = std::max(prevp, satdt / static_cast<Real>(2.0));
                } else {
                    prevp = std::min(prevp, satdt / static_cast<Real>(2.0));
                }
            }

            Real supcol = t0c - t;
            Real n0sfac_val = std::max(std::min(libm::exp(alpha * supcol), n0smax / n0s), static_cast<Real>(1.0));
            Real supsat_ice = std::max(q, qmin) - qsat[1];
            Real satdt_ice = supsat_ice / dtcld;
            int ifsat = 0;

            Real temp = den * std::max(qi, qmin);
            temp = std::sqrt(std::sqrt(temp * temp * temp));
            xni = std::min(std::max(static_cast<Real>(5.38e7) * temp, static_cast<Real>(1.0e3)), static_cast<Real>(1.0e6));
            Real eacrs = libm::exp(wlit<Real>(0.07) * (-supcol));

            Real xmi = den * qi / xni;
            Real diameter = std::min(dicon * std::sqrt(xmi), dimax);
            Real vt2i = wlit<Real>(1.49e4) * libm::pow(diameter, static_cast<Real>(1.31));
            Real vt2r = pvtr * rslopeb[0] * denfac;
            Real vt2s = pvts * rslopeb[1] * denfac;
            Real vt2g = pvtg * rslopeb[2] * denfac;

            Real qsum_val = std::max(qs + qg, static_cast<Real>(1.0e-15));
            Real vt2ave = (qsum_val > static_cast<Real>(1.0e-15))
                ? (vt2s * qs + vt2g * qg) / qsum_val : static_cast<Real>(0.0);

            Real praci = wlit<Real>(0.0);
            Real piacr = wlit<Real>(0.0);
            Real psaci = wlit<Real>(0.0);
            Real pgaci = wlit<Real>(0.0);

            if (supcol > wlit<Real>(0.0) && qi > qmin) {
                if (qr > qcrmin) {
                    Real acrfac = wlit<Real>(2.0) * rslope3[0] + wlit<Real>(2.0) * diameter * rslope2[0]
                        + diameter * diameter * rslope[0];
                    praci = pi * qi * n0r * std::abs(vt2r - vt2i) * acrfac / wlit<Real>(4.0);
                    Real eff = std::min(std::max(static_cast<Real>(0.0), qr / qi), static_cast<Real>(1.0));
                    praci = praci * (eff * eff);
                    praci = std::min(praci, qi / dtcld);

                    piacr = (pi * pi) * avtr * n0r * denr * xni * denfac
                        * g6pbr * rslope3[0] * rslope3[0]
                        * rslopeb[0] / wlit<Real>(24.0) / den;
                    Real eff2 = std::min(std::max(static_cast<Real>(0.0), qi / qr), static_cast<Real>(1.0));
                    piacr = piacr * (eff2 * eff2);
                    piacr = std::min(piacr, qr / dtcld);
                }
                if (qs > qcrmin) {
                    Real acrfac = wlit<Real>(2.0) * rslope3[1] + wlit<Real>(2.0) * diameter * rslope2[1]
                        + diameter * diameter * rslope[1];
                    psaci = pi * qi * eacrs * n0s * n0sfac_val * std::abs(vt2ave - vt2i) * acrfac / wlit<Real>(4.0);
                    psaci = std::min(psaci, qi / dtcld);
                }
                if (qg > qcrmin) {
                    Real egi = libm::exp(wlit<Real>(0.07) * (-supcol));
                    Real acrfac = wlit<Real>(2.0) * rslope3[2] + wlit<Real>(2.0) * diameter * rslope2[2]
                        + diameter * diameter * rslope[2];
                    pgaci = pi * egi * qi * n0g * std::abs(vt2ave - vt2i) * acrfac / wlit<Real>(4.0);
                    pgaci = std::min(pgaci, qi / dtcld);
                }
            }

            Real psacw = wlit<Real>(0.0);
            if (qs > qcrmin && qc > qmin) {
                Real eff = std::min(std::max(static_cast<Real>(0.0), qs / qc), static_cast<Real>(1.0));
                psacw = std::min(pacrc * n0sfac_val * rslope3[1] * rslopeb[1]
                    * (eff * eff) * qc * denfac, qc / dtcld);
            }

            Real pgacw = wlit<Real>(0.0);
            if (qg > qcrmin && qc > qmin) {
                Real eff = std::min(std::max(static_cast<Real>(0.0), qg / qc), static_cast<Real>(1.0));
                pgacw = std::min(pacrg * rslope3[2] * rslopeb[2]
                    * (eff * eff) * qc * denfac, qc / dtcld);
            }

            Real paacw = wlit<Real>(0.0);
            if (qsum_val > static_cast<Real>(1.0e-15)) {
                paacw = (qs * psacw + qg * pgacw) / qsum_val;
            }

            Real pracs = wlit<Real>(0.0);
            Real psacr = wlit<Real>(0.0);
            if (qs > qcrmin && qr > qcrmin) {
                if (supcol > wlit<Real>(0.0)) {
                    Real acrfac = wlit<Real>(5.0) * rslope3[1] * rslope3[1] * rslope[0]
                        + wlit<Real>(2.0) * rslope3[1] * rslope2[1] * rslope2[0]
                        + wlit<Real>(0.5) * rslope2[1] * rslope2[1] * rslope3[0];
                    pracs = (pi * pi) * n0r * n0s * n0sfac_val * std::abs(vt2r - vt2ave)
                        * (dens / den) * acrfac;
                    Real eff = std::min(std::max(static_cast<Real>(0.0), qr / qs), static_cast<Real>(1.0));
                    pracs = pracs * (eff * eff);
                    pracs = std::min(pracs, qs / dtcld);
                }
                Real acrfac = wlit<Real>(5.0) * rslope3[0] * rslope3[0] * rslope[1]
                    + wlit<Real>(2.0) * rslope3[0] * rslope2[0] * rslope2[1]
                    + wlit<Real>(0.5) * rslope2[0] * rslope2[0] * rslope3[1];
                psacr = (pi * pi) * n0r * n0s * n0sfac_val * std::abs(vt2ave - vt2r)
                    * (denr / den) * acrfac;
                Real eff = std::min(std::max(static_cast<Real>(0.0), qs / qr), static_cast<Real>(1.0));
                psacr = psacr * (eff * eff);
                psacr = std::min(psacr, qr / dtcld);
            }

            Real pgacr = wlit<Real>(0.0);
            if (qg > qcrmin && qr > qcrmin) {
                Real acrfac = wlit<Real>(5.0) * rslope3[0] * rslope3[0] * rslope[2]
                    + wlit<Real>(2.0) * rslope3[0] * rslope2[0] * rslope2[2]
                    + wlit<Real>(0.5) * rslope2[0] * rslope2[0] * rslope3[2];
                pgacr = (pi * pi) * n0r * n0g * std::abs(vt2ave - vt2r) * (denr / den) * acrfac;
                Real eff = std::min(std::max(static_cast<Real>(0.0), qg / qr), static_cast<Real>(1.0));
                pgacr = pgacr * (eff * eff);
                pgacr = std::min(pgacr, qr / dtcld);
            }

            Real pgacs = wlit<Real>(0.0);
            Real pseml = wlit<Real>(0.0);
            Real pgeml = wlit<Real>(0.0);
            if (supcol <= wlit<Real>(0.0)) {
                Real xlf = xlf0;
                if (qs > wlit<Real>(0.0)) {
                    pseml = std::min(std::max(cliq * supcol * (paacw + psacr) / xlf, -qs / dtcld), static_cast<Real>(0.0));
                }
                if (qg > wlit<Real>(0.0)) {
                    pgeml = std::min(std::max(cliq * supcol * (paacw + pgacr) / xlf, -qg / dtcld), static_cast<Real>(0.0));
                }
            }

            Real pidep = wlit<Real>(0.0);
            Real psdep = wlit<Real>(0.0);
            Real pgdep = wlit<Real>(0.0);
            Real pigen = wlit<Real>(0.0);
            Real psaut = wlit<Real>(0.0);
            Real pgaut = wlit<Real>(0.0);

            if (supcol > wlit<Real>(0.0)) {
                if (qi > wlit<Real>(0.0) && ifsat != 1) {
                    pidep = wlit<Real>(4.0) * diameter * xni * (rh[1] - wlit<Real>(1.0)) / work1_1;
                    Real supice = satdt_ice - prevp;
                    if (pidep < wlit<Real>(0.0)) {
                        pidep = std::max(std::max(pidep, satdt_ice / static_cast<Real>(2.0)), supice);
                        pidep = std::max(pidep, -qi / dtcld);
                    } else {
                        pidep = std::min(std::min(pidep, satdt_ice / static_cast<Real>(2.0)), supice);
                    }
                    if (std::abs(prevp + pidep) >= std::abs(satdt_ice)) ifsat = 1;
                }

                if (qs > wlit<Real>(0.0) && ifsat != 1) {
                    Real coeres = rslope2[1] * std::sqrt(rslope[1] * rslopeb[1]);
                    psdep = (rh[1] - wlit<Real>(1.0)) * n0sfac_val * (precs1 * rslope2[1]
                        + precs2 * work2_k * coeres) / work1_1;
                    Real supice = satdt_ice - prevp - pidep;
                    if (psdep < wlit<Real>(0.0)) {
                        psdep = std::max(psdep, -qs / dtcld);
                        psdep = std::max(std::max(psdep, satdt_ice / static_cast<Real>(2.0)), supice);
                    } else {
                        psdep = std::min(std::min(psdep, satdt_ice / static_cast<Real>(2.0)), supice);
                    }
                    if (std::abs(prevp + pidep + psdep) >= std::abs(satdt_ice)) ifsat = 1;
                }

                if (qg > wlit<Real>(0.0) && ifsat != 1) {
                    Real coeres = rslope2[2] * std::sqrt(rslope[2] * rslopeb[2]);
                    pgdep = (rh[1] - wlit<Real>(1.0)) * (precg1 * rslope2[2]
                        + precg2 * work2_k * coeres) / work1_1;
                    Real supice = satdt_ice - prevp - pidep - psdep;
                    if (pgdep < wlit<Real>(0.0)) {
                        pgdep = std::max(pgdep, -qg / dtcld);
                        pgdep = std::max(std::max(pgdep, satdt_ice / static_cast<Real>(2.0)), supice);
                    } else {
                        pgdep = std::min(std::min(pgdep, satdt_ice / static_cast<Real>(2.0)), supice);
                    }
                    if (std::abs(prevp + pidep + psdep + pgdep) >= std::abs(satdt_ice)) ifsat = 1;
                }

                if (supsat_ice > wlit<Real>(0.0) && ifsat != 1) {
                    Real supice = satdt_ice - prevp - pidep - psdep - pgdep;
                    Real xni0 = wlit<Real>(1.0e3) * libm::exp(wlit<Real>(0.1) * supcol);
                    Real roqi0 = wlit<Real>(4.92e-11) * libm::pow(xni0, static_cast<Real>(1.33));
                    pigen = std::max(static_cast<Real>(0.0), (roqi0 / den - std::max(qi, static_cast<Real>(0.0))) / dtcld);
                    pigen = std::min(std::min(pigen, satdt_ice), supice);
                }

                if (qi > wlit<Real>(0.0)) {
                    Real qimax = roqimax / den;
                    psaut = std::max(static_cast<Real>(0.0), (qi - qimax) / dtcld);
                }

                if (qs > wlit<Real>(0.0)) {
                    Real alpha2 = wlit<Real>(1.0e-3) * libm::exp(wlit<Real>(0.09) * (-supcol));
                    pgaut = std::min(std::max(static_cast<Real>(0.0), alpha2 * (qs - qs0)), qs / dtcld);
                }
            }

            Real psevp = wlit<Real>(0.0);
            Real pgevp = wlit<Real>(0.0);
            if (supcol < wlit<Real>(0.0)) {
                if (qs > wlit<Real>(0.0) && rh[0] < wlit<Real>(1.0)) {
                    Real coeres = rslope2[1] * std::sqrt(rslope[1] * rslopeb[1]);
                    psevp = (rh[0] - wlit<Real>(1.0)) * n0sfac_val * (precs1 * rslope2[1]
                        + precs2 * work2_k * coeres) / work1_0;
                    psevp = std::min(std::max(psevp, -qs / dtcld), static_cast<Real>(0.0));
                }
                if (qg > wlit<Real>(0.0) && rh[0] < wlit<Real>(1.0)) {
                    Real coeres = rslope2[2] * std::sqrt(rslope[2] * rslopeb[2]);
                    pgevp = (rh[0] - wlit<Real>(1.0)) * (precg1 * rslope2[2]
                        + precg2 * work2_k * coeres) / work1_0;
                    pgevp = std::min(std::max(pgevp, -qg / dtcld), static_cast<Real>(0.0));
                }
            }

            // Mass conservation check and feedback
            Real delta2 = (qr < wlit<Real>(1.0e-4) && qs < wlit<Real>(1.0e-4)) ? static_cast<Real>(1.0) : static_cast<Real>(0.0);
            Real delta3 = (qr < wlit<Real>(1.0e-4)) ? static_cast<Real>(1.0) : static_cast<Real>(0.0);

            if (t <= t0c) {
                // cloud water
                Real value = std::max(qmin, qc);
                Real source = (praut + pracw + paacw + paacw) * dtcld;
                if (source > value) {
                    Real factor = value / source;
                    praut *= factor;
                    pracw *= factor;
                    paacw *= factor;
                }

                // cloud ice
                value = std::max(qmin, qi);
                source = (psaut - pigen - pidep + praci + psaci + pgaci) * dtcld;
                if (source > value) {
                    Real factor = value / source;
                    psaut *= factor;
                    pigen *= factor;
                    pidep *= factor;
                    praci *= factor;
                    psaci *= factor;
                    pgaci *= factor;
                }

                // rain
                value = std::max(qmin, qr);
                source = (-praut - prevp - pracw + piacr + psacr + pgacr) * dtcld;
                if (source > value) {
                    Real factor = value / source;
                    praut *= factor;
                    prevp *= factor;
                    pracw *= factor;
                    piacr *= factor;
                    psacr *= factor;
                    pgacr *= factor;
                }

                // snow
                value = std::max(qmin, qs);
                source = -(psdep + psaut - pgaut + paacw + piacr * delta3
                    + praci * delta3 - pracs * (wlit<Real>(1.0) - delta2)
                    + psacr * delta2 + psaci - pgacs) * dtcld;
                if (source > value) {
                    Real factor = value / source;
                    psdep *= factor;
                    psaut *= factor;
                    pgaut *= factor;
                    paacw *= factor;
                    piacr *= factor;
                    praci *= factor;
                    psaci *= factor;
                    pracs *= factor;
                    psacr *= factor;
                    pgacs *= factor;
                }

                // graupel
                value = std::max(qmin, qg);
                source = -(pgdep + pgaut + piacr * (wlit<Real>(1.0) - delta3) + praci * (wlit<Real>(1.0) - delta3)
                    + psacr * (wlit<Real>(1.0) - delta2) + pracs * (wlit<Real>(1.0) - delta2)
                    + pgaci + paacw + pgacr + pgacs) * dtcld;
                if (source > value) {
                    Real factor = value / source;
                    pgdep *= factor;
                    pgaut *= factor;
                    piacr *= factor;
                    praci *= factor;
                    psacr *= factor;
                    pracs *= factor;
                    paacw *= factor;
                    pgaci *= factor;
                    pgacr *= factor;
                    pgacs *= factor;
                }

                Real work2_val = -(prevp + psdep + pgdep + pigen + pidep);
                q += work2_val * dtcld;
                qc = std::max(qc - (praut + pracw + paacw + paacw) * dtcld, static_cast<Real>(0.0));
                qr = std::max(qr + (praut + pracw + prevp - piacr - pgacr - psacr) * dtcld, static_cast<Real>(0.0));
                qi = std::max(qi - (psaut + praci + psaci + pgaci - pigen - pidep) * dtcld, static_cast<Real>(0.0));
                qs = std::max(qs + (psdep + psaut + paacw - pgaut + piacr * delta3
                    + praci * delta3 + psaci - pgacs - pracs * (static_cast<Real>(1.0) - delta2) + psacr * delta2) * dtcld, static_cast<Real>(0.0));
                qg = std::max(qg + (pgdep + pgaut + piacr * (static_cast<Real>(1.0) - delta3)
                    + praci * (static_cast<Real>(1.0) - delta3) + psacr * (static_cast<Real>(1.0) - delta2)
                    + pracs * (static_cast<Real>(1.0) - delta2) + pgaci + paacw + pgacr + pgacs) * dtcld, static_cast<Real>(0.0));

                Real xlf = xls - xl;
                Real xlwork2 = -xls * (psdep + pgdep + pidep + pigen)
                    - xl * prevp - xlf * (piacr + paacw + paacw + pgacr + psacr);
                t -= xlwork2 / cpm * dtcld;

            } else {
                // t > t0c
                // cloud water
                Real value = std::max(qmin, qc);
                Real source = (praut + pracw + paacw + paacw) * dtcld;
                if (source > value) {
                    Real factor = value / source;
                    praut *= factor;
                    pracw *= factor;
                    paacw *= factor;
                }

                // rain
                value = std::max(qmin, qr);
                source = (-paacw - praut + pseml + pgeml - pracw - paacw - prevp) * dtcld;
                if (source > value) {
                    Real factor = value / source;
                    praut *= factor;
                    prevp *= factor;
                    pracw *= factor;
                    paacw *= factor;
                    pseml *= factor;
                    pgeml *= factor;
                }

                // snow
                value = std::max(qcrmin, qs);
                source = (pgacs - pseml - psevp) * dtcld;
                if (source > value) {
                    Real factor = value / source;
                    pgacs *= factor;
                    psevp *= factor;
                    pseml *= factor;
                }

                // graupel
                value = std::max(qcrmin, qg);
                source = -(pgacs + pgevp + pgeml) * dtcld;
                if (source > value) {
                    Real factor = value / source;
                    pgacs *= factor;
                    pgevp *= factor;
                    pgeml *= factor;
                }

                Real work2_val = -(prevp + psevp + pgevp);
                q += work2_val * dtcld;
                qc = std::max(qc - (praut + pracw + paacw + paacw) * dtcld, static_cast<Real>(0.0));
                qr = std::max(qr + (praut + pracw + prevp + paacw + paacw - pseml - pgeml) * dtcld, static_cast<Real>(0.0));
                qs = std::max(qs + (psevp - pgacs + pseml) * dtcld, static_cast<Real>(0.0));
                qg = std::max(qg + (pgacs + pgevp + pgeml) * dtcld, static_cast<Real>(0.0));

                Real xlf = xls - xl;
                Real xlwork2 = -xl * (prevp + psevp + pgevp) - xlf * (pseml + pgeml);
                t -= xlwork2 / cpm * dtcld;
            }
        }
        Real hsub = xls;
        Real hvap = xlv0;
        Real cvap = cpv;
        Real ttp = t0c + wlit<Real>(0.01);
        Real dldt = cvap - cliq;
        Real xa = -dldt / rv;
        Real xb = xa + hvap / (rv * ttp);
        Real dldti = cvap - cice;
        Real xai = -dldti / rv;
        Real xbi = xai + hsub / (rv * ttp);
        (void)xbi;
        {
            Real tr = ttp / t;
            Real qs1 = psat * libm::exp(libm::log(tr) * xa) * libm::exp(xb * (wlit<Real>(1.0) - tr));
            qs1 = std::min(qs1, static_cast<Real>(0.99) * p);
            qs1 = ep2 * qs1 / (p - qs1);
            qs1 = std::max(qs1, qmin);

            Real work1_cond = conden(t, q, qs1, xl, cpm);
            Real pcond = std::min(std::max(work1_cond / dtcld, static_cast<Real>(0.0)),
                                  std::max(q, static_cast<Real>(0.0)) / dtcld);
            if (qc > wlit<Real>(0.0) && work1_cond < wlit<Real>(0.0)) {
                pcond = std::max(work1_cond, -qc) / dtcld;
            }
            q -= pcond * dtcld;
            qc = std::max(qc + pcond * dtcld, static_cast<Real>(0.0));
            t += pcond * xl / cpm * dtcld;
        }
        if (qc <= qmin) qc = wlit<Real>(0.0);
        if (qi <= qmin) qi = wlit<Real>(0.0);
    }
#endif
};

}  // namespace cppmpas::phys::wsm6
