// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// src/operators/mpas_geometry_utils.F: the functions atm_initialize_deformation_weights calls.
#include "core/precision.hpp"

namespace cppmpas {

// mpas_geometry_utils.F:131-154
real mpas_arc_length(real ax, real ay, real az, real bx, real by, real bz);
// mpas_geometry_utils.F:27-72
real mpas_sphere_angle(real ax, real ay, real az, real bx, real by, real bz, real cx, real cy, real cz);
// mpas_geometry_utils.F:81-121
real mpas_plane_angle(real ax, real ay, real az, real bx, real by, real bz, real cx, real cy, real cz, real u, real v,
                      real w);

// mpas_vector_operations.F:861-889
real mpas_fix_periodicity(real pxi, real xci, real xiRef);

// gfortran's MAX(a, b) / MIN(a, b) as the reference executable computes them here: maxsd / minsd with a in the
// destination, i.e. b unless a > b (a < b); the result for NaN and for +-0 ties is b.
inline real maxsd(real a, real b) { return a > b ? a : b; }
inline real minsd(real a, real b) { return a < b ? a : b; }

}  // namespace cppmpas
