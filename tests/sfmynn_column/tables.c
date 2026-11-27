// SPDX-License-Identifier: GPL-3.0-or-later
// Copies sf_mynn's private module tables (MMM-physics sf_mynn.F90:81-82, the global symbols gfortran emits for them)
// for harness.F90; no arithmetic.
#include <string.h>
extern double __sf_mynn_MOD_psim_stab[1001], __sf_mynn_MOD_psim_unstab[1001], __sf_mynn_MOD_psih_stab[1001],
    __sf_mynn_MOD_psih_unstab[1001];
void get_sfmynn_tables(double* t) {
  memcpy(t, __sf_mynn_MOD_psim_stab, sizeof(double) * 1001);
  memcpy(t + 1001, __sf_mynn_MOD_psim_unstab, sizeof(double) * 1001);
  memcpy(t + 2002, __sf_mynn_MOD_psih_stab, sizeof(double) * 1001);
  memcpy(t + 3003, __sf_mynn_MOD_psih_unstab, sizeof(double) * 1001);
}
