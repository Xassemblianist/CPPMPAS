// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// thompson_init (physics/physics_wrf/module_mp_thompson.F, MT:433-993) as MPAS's init_microphysics calls it
// (mpas_atmphys_driver_microphysics.F:295, l_mp_tables = .true., mpas_atmphys_control.F:406): the module constants,
// the size bins, the collision-efficiency tables (table_Efrw MT:4068-4124, table_Efsw MT:4132-4168), the
// cloud-drop evaporation table (table_dropEvap MT:4222-4310), the four lookup-table files MPAS's build_tables
// writes (MP_THOMPSON_QRacrQG / QRacrQS / freezeH2O / QIautQS_DATA.DBL, MT:850-965), and radar_init with
// Thompson's mass / shape parameters (MT:978-988).
//
// The tables are read from the reference run's own .DBL files (big-endian Fortran unformatted sequential records,
// -fconvert=big-endian, as the reference's thompson_init reads them); build_tables
// (utils/atmphys_build_tables_thompson.F) is not ported. tests/thompson_column checks every array and constant
// thompson_init leaves (and the tables read) against MPAS's own thompson_init built with the reference flags.
//
// Reference object code (__module_mp_thompson_MOD_thompson_init, 4635 instructions; the standalone -O3 build is
// instruction-identical): 38 WGAMMA calls (mpas_atmphys_functions, scalar log / exp), 1 log (nic1), 10
// _ZGVbN2v_exp (the bin loops MT:684, 696, 708, 720, 732, 2 per loop), no pow: gfortran folds the powers of
// constants with MPFR (Sc3, D0i, xm0g, ocms, ocmg), the logs of the bin loops, SQRT(av_x), and NINT(ALOG10(..)).
// table_Efrw / table_Efsw: scalar log, exp, atan, pow; table_dropEvap: _ZGVbN2vv_pow for massc (MT:4235, 100
// values, 50 pairs), scalar pow / exp, __powidf2 for Dc**nu_c.
#include <cstdint>
#include <string>
#include <vector>

#include "core/field.hpp"
#include "physics/mp_thompson.hpp"

namespace cppmpas::phys::thompson {

// The tables of thompson_init and the view the kernels read (pointers into the vectors below).
struct ThompsonTables {
  ThompsonView v{};
  FStorage<double> tcg_racg, tmr_racg, tcr_gacr, tmg_gacr, tnr_racg, tnr_gacr;
  FStorage<double> tcs_racs1, tmr_racs1, tcs_racs2, tmr_racs2, tcr_sacr1, tms_sacr1, tcr_sacr2, tms_sacr2,
      tnr_racs1, tnr_racs2, tnr_sacr1, tnr_sacr2;
  FStorage<double> tpi_qrfz, tni_qrfz, tpg_qrfz, tnr_qrfz, tpi_qcfz, tni_qcfz;
  FStorage<double> tpi_ide, tps_iaus, tni_iaus;
  FStorage<double> t_Efrw, t_Efsw;
  FStorage<double> tpc_wev, tnc_wev;  // (nbc, ntb_c, nbc): read only by mp_thompson_aerosols
  double xDx[nbins + 1];                 // the last bin edges (cloud droplet number), a module array
  void point();                          // sets the view's table pointers to the vectors
};

// thompson_init(l_mp_tables = .true.); dir: the directory of the four MP_THOMPSON_*_DATA.DBL files. Errors as
// MPAS's physics_error_fatal ('failure opening ...') and a short or malformed record.
void thompson_init(const std::string& dir, ThompsonTables& t);

// Every array and constant thompson_init leaves, in a fixed order (name, values), for the column check's
// comparison with MPAS's own thompson_init. Tables over 100000 values are reduced to their FNV-1a hash.
void dump_thompson_init(const ThompsonTables& t, const std::string& path);

}  // namespace cppmpas::phys::thompson
