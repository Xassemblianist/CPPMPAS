// SPDX-License-Identifier: GPL-3.0-or-later
#include "operators/geometry_utils.hpp"

#include <cmath>

#include "core/libm.hpp"

namespace cppmpas {

// mpas_geometry_utils.F:131-154. In the reference (inlined into mpas_sphere_angle as well) 2.0*r is r+r and
// r*2.0*asin(...) is asin(...)*(r+r): both exact rewrites.
real mpas_arc_length(real ax, real ay, real az, real bx, real by, real bz) {
  const real cx = bx - ax;  // :140-142
  const real cy = by - ay;
  const real cz = bz - az;
  const real r = std::sqrt(ax * ax + ay * ay + az * az);  // :149
  const real c = std::sqrt(cx * cx + cy * cy + cz * cz);  // :150
  return r * 2.0 * libm::asin(c / (2.0 * r));             // :152, asin@plt (__ieee754_asin_fma)
}

// mpas_geometry_utils.F:27-72
real mpas_sphere_angle(real ax, real ay, real az, real bx, real by, real bz, real cx, real cy, real cz) {
  const real a = mpas_arc_length(bx, by, bz, cx, cy, cz);  // :46-48
  const real b = mpas_arc_length(ax, ay, az, cx, cy, cz);
  const real c = mpas_arc_length(ax, ay, az, bx, by, bz);

  const real ABx = bx - ax;  // :50-56
  const real ABy = by - ay;
  const real ABz = bz - az;
  const real ACx = cx - ax;
  const real ACy = cy - ay;
  const real ACz = cz - az;

  const real Dx = (ABy * ACz) - (ABz * ACy);     // :58-60
  const real Dy = -((ABx * ACz) - (ABz * ACx));
  const real Dz = (ABx * ACy) - (ABy * ACx);

  const real s = 0.5 * (a + b + c);  // :62
  // :64, sin@plt x4 (__sin_fma); MAX / MIN as maxsd / minsd (operand order of the reference executable)
  const real q = (libm::sin(s - b) * libm::sin(s - c)) / (libm::sin(b) * libm::sin(c));
  const real sin_angle = std::sqrt(minsd(maxsd(q, 0.0), 1.0));

  if ((Dx * ax + Dy * ay + Dz * az) >= 0.0) {  // :66-70
    return 2.0 * libm::asin(maxsd(minsd(sin_angle, 1.0), -1.0));
  } else {
    return -2.0 * libm::asin(maxsd(minsd(sin_angle, 1.0), -1.0));
  }
}

// mpas_geometry_utils.F:81-121. x**2.0 is x*x (gcc folds pow(x, 2.0)). In atm_initialize_deformation_weights the
// result is overwritten before use (mpas_atm_core.F:1807 vs 1821), but the call is made, as in the reference.
real mpas_plane_angle(real ax, real ay, real az, real bx, real by, real bz, real cx, real cy, real cz, real u, real v,
                      real w) {
  const real ABx = bx - ax;  // :98-101
  const real ABy = by - ay;
  const real ABz = bz - az;
  const real mAB = std::sqrt(ABx * ABx + ABy * ABy + ABz * ABz);
  const real ACx = cx - ax;  // :103-106
  const real ACy = cy - ay;
  const real ACz = cz - az;
  const real mAC = std::sqrt(ACx * ACx + ACy * ACy + ACz * ACz);
  const real Dx = (ABy * ACz) - (ABz * ACy);  // :109-111
  const real Dy = -((ABx * ACz) - (ABz * ACx));
  const real Dz = (ABx * ACy) - (ABy * ACx);
  const real cos_angle = (ABx * ACx + ABy * ACy + ABz * ACz) / (mAB * mAC);  // :113
  if ((Dx * u + Dy * v + Dz * w) >= 0.0) {                                   // :115-119, acos@plt
    return libm::acos(maxsd(minsd(cos_angle, 1.0), -1.0));
  } else {
    return -libm::acos(maxsd(minsd(cos_angle, 1.0), -1.0));
  }
}

// mpas_vector_operations.F:861-889
real mpas_fix_periodicity(real pxi, real xci, real xiRef) {
  const real dist = pxi - xci;
  if (std::fabs(dist) > xiRef * 0.5) return pxi - (dist / std::fabs(dist)) * xiRef;
  return pxi;
}

}  // namespace cppmpas
