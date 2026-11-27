// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// GAMMLN and WGAMMA of physics/mpas_atmphys_functions.F:98-120, 144-151 (MPAS v8.4.2), the gamma functions
// Thompson's thompson_init calls (module_mp_thompson.F uses mpas_atmphys_functions' gammln / wgamma with -Dmpas).
// Reference object code (__mpas_atmphys_functions_MOD_gammln / _wgamma; the standalone -O3 build is
// instruction-identical): gammln 2 log@plt, wgamma = EXP(GAMMLN(y)) with gammln inlined, 2 log@plt + exp@plt; no FMA.
// RSLF / RSIF of the same file are in physics/mpas_atmphys_microphysics.cpp.
#include "core/backend.hpp"
#include "core/libm.hpp"

namespace cppmpas::phys {

// mpas_atmphys_functions.F:98-120
CPPMPAS_HD inline double gammln(double xx) {
  constexpr double STP = 2.5066282746310005;
  constexpr double COF[6] = {76.18009172947146, -86.50532032941677, 24.01409824083091,
                             -1.231739572450155, .1208650973866179e-2, -.5395239384953e-5};
  const double X = xx;
  double Y = X;
  double TMP = X + 5.5;
  TMP = (X + 0.5) * ::cppmpas::libm::log(TMP) - TMP;
  double SER = 1.000000000190015;
  for (int J = 0; J < 6; ++J) {
    Y = Y + 1.;
    SER = SER + COF[J] / Y;
  }
  return TMP + ::cppmpas::libm::log(STP * SER / X);
}

// mpas_atmphys_functions.F:144-151
CPPMPAS_HD inline double wgamma(double y) { return ::cppmpas::libm::exp(gammln(y)); }

}  // namespace cppmpas::phys
