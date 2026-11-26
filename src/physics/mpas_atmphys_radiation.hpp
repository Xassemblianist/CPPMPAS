// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// MPAS-Atmosphere's radiation drivers for mesoscale_reference (config_radt_lw_scheme = 'rrtmg_lw',
// config_radt_sw_scheme = 'rrtmg_sw', config_radt_cld_scheme = 'cld_fraction', mpas_atmphys_control.F:139-141),
// statement by statement, on the arrays of the MPAS pools (Fortran layout, 1-based) and the module arrays of
// mpas_atmphys_vars. One rank, one thread: its..ite = 1..nCellsSolve, kts..kte = 1..nVertLevels, jts = jte = 1.
//  * physics/mpas_atmphys_driver_cloudiness.F (CL): allocate / deallocate_cloudiness (CL:64-89),
//    cloudiness_from_MPAS (CL:92-143), cloudiness_to_MPAS (CL:146-172), driver_cloudiness (CL:175-234) with
//    calc_cldfraction (CL:284-370);
//  * physics/mpas_atmphys_driver_radiation_sw.F (RS): radiation_sw_from_MPAS (RS:355-711), radiation_sw_to_MPAS
//    (RS:714-794), driver_radiation_sw (RS:834-1006), radconst (RS:1009-1050), init_radiation_sw (RS:797-831);
//  * physics/mpas_atmphys_driver_radiation_lw.F (RL): radiation_lw_from_MPAS (RL:297-610), radiation_lw_to_MPAS
//    (RL:613-718), driver_radiation_lw (RL:804-960), init_radiation_lw (RL:768-801);
//  * physics/mpas_atmphys_o3climatology.F (OZ): init_o3climatology (OZ:62-180), o3climatology_from_MPAS
//    (OZ:228-325) with module_ra_cam_support.F's getfactors (physics_wrf/module_ra_cam_support.F:2528-2611);
//  * physics/mpas_atmphys_manager.F (MG): the radiation time variables of physics_run_init (MG:457-461) and
//    physics_timetracker (MG:199-206).
// Refused (hard errors): 'cam_lw' / 'cam_sw' and any other scheme, 'cld_incidence' and 'cld_fraction_thompson',
// config_microp_re = true (has_reqc / has_reqi / has_reqs), 'mp_thompson_aerosols' (aer_opt = 3),
// config_radt_cld_overlap = 'exponential' (a libmvec lane, module_ra_rrtmg.hpp).
//
// Reference object code (bin/atmosphere_model, sha256 524a1cf4...): no FMA, no ymm register. radconst calls asin,
// sin and sincos (cos(rjul) / sin(rjul) and cos(2rjul) / sin(2rjul) as two sincos); driver_cloudiness (with
// calc_cldfraction and cloudiness_from_MPAS inlined) calls exp three times and pow, scalar, and one _ZGVbN2vv_pow:
// the dx_p loop of CL:122-127 (len_disp / meshDensity**0.25), vectorised over i. dx_p is used only by
// 'cld_fraction_thompson' and deallocate_cloudiness frees it before any other scheme reads it (the surface layer
// computes its own), so that loop is left out.
#include <string>
#include <vector>

#include "core/field.hpp"
#include "core/precision.hpp"
#include "physics/module_ra_rrtmg.hpp"
#include "physics/mpas_atmphys_interface.hpp"

namespace cppmpas::phys {

struct RadiationPools {
  int nVertLevels = 0;
  int nOznLevels = 0, nMonths = 0;  // mesh dimensions (num_oznlevels, num_months)
  // configs
  std::string radt_lw_scheme, radt_sw_scheme, radt_cld_scheme, radt_cld_overlap, radt_cld_dcorrlen, microp_scheme;
  bool config_o3climatology = true, config_microp_re = false;
  // mesh
  FView<const real> latCell, lonCell;
  // sfc_input
  FView<const real> skintemp, snow, xice, xland;
  // atm_input (filled by init_o3climatology)
  FView<real> pin, ozmixm;  // pin(nOznLevels), ozmixm(nMonths, nOznLevels, nCells)
  // diag_physics
  FView<const real> sfc_albedo, sfc_emiss;
  FView<real> cldfrac, o3clim, o3vmr, rre_cloud, rre_ice, rre_snow;
  FView<real> glw, lwcf, lwdnb, lwdnbc, lwdnt, lwdntc, lwupb, lwupbc, lwupt, lwuptc, olrtoa;
  FView<real> coszr, gsw, swcf, swdnb, swdnbc, swdnt, swdntc, swupb, swupbc, swupt, swuptc, swddir, swddni, swddif;
  // tend_physics
  FView<real> rthratenlw, rthratensw;
};

// The radiation variables of mpas_atmphys_vars and mpas_atmphys_manager, and the module arrays the cloudiness and
// radiation drivers share within a physics step (allocated by allocate_cloudiness / allocate_radiation_sw / _lw,
// freed after the step: mpas_atmphys_driver.F:202-246)
struct RadiationVars {
  int year = 0, julday = 0;
  real gmt = 0, curr_julday = 0;       // MG:457-461 (gmt from the start time only), MG:203-206
  real dt_radtlw = 0, dt_radtsw = 0;   // MG:466-490 (config_radt*_interval in seconds)
  real declin = 0, solcon = 0;         // radconst
  int has_reqc = 0, has_reqi = 0, has_reqs = 0;  // MG:747-760 (config_microp_re = false)
  int icloud = 1;                      // mpas_atmphys_vars.F:807
  bool sw_ran = false;                 // xlat_p set by this step's radiation_sw_from_MPAS
  // cloudiness (CL:64-75)
  PArray cldfrac_p, qvrad_p, qcrad_p, qirad_p, qsrad_p;
  std::vector<real> xland_cl;  // CL's xland_p (the same module array as the radiation drivers')
  // radiation (RS:106-235, RL:106-187)
  PArray f_ice, f_rain, rthratenlw_p, rthratensw_p, recloud_p, reice_p, resnow_p, rrecloud_p, rreice_p, rresnow_p,
      o3clim_p;
  std::vector<real> pin_p, xlat_p, xlon_p, sfc_albedo_p, sfc_emiss_p, snow_p, tsk_p, xice_p, xland_p;
  std::vector<real> glw_p, lwcf_p, lwdnb_p, lwdnbc_p, lwdnt_p, lwdntc_p, lwupb_p, lwupbc_p, lwupt_p, lwuptc_p,
      olrtoa_p;
  std::vector<real> coszr_p, gsw_p, swcf_p, swdnb_p, swdnbc_p, swdnt_p, swdntc_p, swupb_p, swupbc_p, swupt_p,
      swuptc_p, swddir_p, swddni_p, swddif_p;
};

// RRTMG tables after rrtmg_initlw_forMPAS / rrtmg_initsw_forMPAS (mpas_atmphys_rrtmg_lwinit.F:43-58,
// _swinit.F): memory images of the tables in CPPWRF's layout (physics/rrtmg/rrtmg_*_tables.hpp)
struct RrtmgTables {
  rrtmg::RrtmgLwTables lw;
  rrtmg::RrtmgSwTables sw;
};
void init_radiation_rrtmg(const std::string& lw_tables_bin, const std::string& sw_tables_bin, RrtmgTables& t);

// MG:457-461 (physics_run_init) and MG:199-206 (physics_timetracker): year, day of year, utc hours from the clock
void physics_radiation_time_init(RadiationVars& r, int yr, int doy, int h, int m, int s, long s_n, long s_d);
void physics_radiation_time_update(RadiationVars& r, int yr, int doy, int h, int m, int s, long s_n, long s_d);

// RS:1009-1050
void radconst(real& declin, real& solcon, real julian, real degrad, real dpd);

// OZ:62-180: reads OZONE_PLEV.TBL, OZONE_LAT.TBL, OZONE_DAT.TBL from dir, fills pin and ozmixm on cells 1..nCells
void init_o3climatology(const std::string& dir, RadiationPools& m, int nCells);
// OZ:228-325 on cells 1..nCellsSolve
void o3climatology_from_MPAS(real julian, RadiationPools& m, int nCellsSolve);

// CL:175-234 (with allocate_cloudiness, CL:64-75)
void driver_cloudiness(RadiationPools& m, const PhysicsVars& v, RadiationVars& r, int its, int ite);
// RS:834-1006 (with allocate_radiation_sw); xtime_s: seconds since the start of the run
void driver_radiation_sw(int itimestep, RadiationPools& m, const PhysicsVars& v, RadiationVars& r,
                         const RrtmgTables& t, real xtime_s, int its, int ite);
// RL:804-960 (with allocate_radiation_lw, RL:103-206)
void driver_radiation_lw(real xtime_s, RadiationPools& m, const PhysicsVars& v, RadiationVars& r,
                         const RrtmgTables& t, int its, int ite);

}  // namespace cppmpas::phys
