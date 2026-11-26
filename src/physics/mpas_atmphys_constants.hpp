// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// The constants of MPAS-Atmosphere's physics: src/framework/mpas_constants.F:46-51 (rgas, rv, cp, rvord, gravity)
// and src/core_atmosphere/physics/mpas_atmphys_constants.F:37-84 (R_d => rgas). With -fdefault-real-8 every literal
// and parameter is a double, and gfortran folds the parameter expressions one correctly rounded operation per node,
// as the constexpr expressions below.
#include "core/precision.hpp"

namespace cppmpas::physics_constants {

constexpr real rgas = 287.0;              // mpas_constants.F:47
constexpr real rv = 461.6;                // :48
constexpr real cp = 7.0 * rgas / 2.0;     // :50
constexpr real rvord = rv / rgas;         // :51
constexpr real gravity = 9.80616;         // :46

constexpr real R_d = rgas;                // mpas_atmphys_constants.F:11
constexpr real P0 = 100000.;              // :40
constexpr real R_v = 461.6;               // :42
constexpr real ep_1 = R_v / R_d - 1.;     // :43
constexpr real ep_2 = R_d / R_v;          // :44
constexpr real cpv = 4. * R_v;            // :45
constexpr real rdg = R_d / gravity;       // :46
constexpr real rcp = R_d / cp;            // :47
constexpr real rcv = R_d / (cp - R_d);    // :48
constexpr real rho_a = 1.28;              // :50
constexpr real rho_r = 1000.0;            // :51
constexpr real rho_s = 100.0;             // :52
constexpr real rho_w = 1000.0;            // :53
constexpr real svp1 = 0.6112;             // :55
constexpr real svp2 = 17.67;              // :56
constexpr real svp3 = 29.65;              // :57
constexpr real svpt0 = 273.15;            // :58
constexpr real xlv = 2.50e6;              // :60
constexpr real xlf = 3.50e5;              // :61
constexpr real xls = xlv + xlf;           // :62
constexpr real karman = 0.4;              // :69
constexpr real cliq = 4190.;              // :73
constexpr real cice = 2106.;              // :74
constexpr real epsilon = 1.e-15;          // :75
constexpr real psat = 610.78;             // :76
constexpr real re_qc_bg = 2.49e-6;        // :78
constexpr real re_qi_bg = 4.99e-6;        // :79
constexpr real re_qs_bg = 9.99e-6;        // :80

}  // namespace cppmpas::physics_constants
