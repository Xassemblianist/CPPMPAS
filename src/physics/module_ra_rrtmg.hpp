// SPDX-License-Identifier: GPL-3.0-or-later AND LicenseRef-AER-RRTMG
#pragma once
// MPAS-Atmosphere's RRTMG column wrappers (MPAS v8.4.2 physics/physics_wrf): rrtmg_lwrad (module_ra_rrtmg_lw.F, LW
// 11574-12305) with inirad / o3data (LW 12308-12444), relcalc (LW 12447-12497) and reicalc (LW 12499-12527),
// rrtmg_swrad (module_ra_rrtmg_sw.F, SW 10025-10667) and vinterp_ozn (module_ra_rrtmg_vinterp.F:16-176), statement by
// statement. They differ from WRF 4.8.0's RRTMG_LWRAD / RRTMG_SWRAD (MPAS rewrote them): the number of buffer layers
// above the model top is set per column from its top pressure (LW 11822), the buffer pressures are bounded by 0
// (LW 11963), there are no time-varying trace gases (the 2005 IPCC values of LW 11681-11700 / SW 10137-10148), the
// CAM ozone climatology is interpolated to every layer, buffer layers included (vinterp_ozn, LW 12030, SW 10379),
// and the SW solar constant comes with its eccentricity factor (radconst; dyofyr = 0, adjes = 1, SW 10396-10397).
//
// The RRTMG kernels below the wrappers (inatm, cldprmc, setcoef, taumol, rtrnmc for the LW; inatm_sw, setcoef_sw,
// taumol_sw, spcvmc_sw with reftra_sw / vrtqdr_sw, cldprmc_sw for the SW; mcica_subcol_lw / _sw with kissvec and
// generate_stochastic_clouds) are CPPWRF's (src/physics/rrtmg, SHARED_WITH_CPPWRF.md): MPAS's
// module_ra_rrtmg_lw.F / _sw.F are WRF 4.8.0's phys/ files up to the wrapper modules (diff: only the RWORDSIZE
// preprocessor block, trailing blanks and wrf_error_fatal spelled physics_error_fatal).
//
// One source for the CPU and the GPU (core/backend.hpp): the columns its..ite go through CPPWRF's RRTMG stages as
// CPPWRF's GPU batches run them (CPPWRF 7fac5bf src/phys/ra_rrtmg_lw/rrtmg_lw_batch.cpp, ra_rrtmg_sw/
// rrtmg_sw_batch.cpp): per column the wrapper's statements up to the call of rrtmg_lw / rrtmg_sw (here
// rrtmg_lwrad_column / rrtmg_swrad_column), then setcoef per (column, layer), taumol per (column, layer), the
// radiative transfer per (column, g-point), its sums per (column, level), and per column the heating rates and the
// wrapper's outputs. The CUDA build launches each stage once for a chunk of columns (one thread per column and layer,
// g-point or level, the column index fastest); the CPU build runs chunks of one column, so every stage loops over that
// column's layers, g-points or levels in order. MPAS parameters, not edits of the shared stages: the per-column number
// of layers (nlay below; the stage launches skip a column's layers above it), MPAS's wrapper statements (instead of
// CPPWRF's rrtmg_lw_prep / rrtmg_sw_prep, which are WRF's) and MPAS's outputs.
//
// Reference object code (bin/atmosphere_model, sha256 524a1cf4...): no FMA and no ymm register in any radiation
// routine. rrtmg_lwrad's only libm call is lround (NINT, LW 11822); rrtmg_swrad calls sincos twice and cos once
// (SW 10205-10206, the zenith angle; xxlat and declin take sincos) . In the kernels every libm call is scalar (exp,
// log, pow) except one libmvec lane: generate_stochastic_clouds' exponential-overlap alpha (_ZGVbN2v_exp, LW 2552 /
// SW 1860), reached only with config_radt_cld_overlap = 'exponential', which is refused here.
#include "core/backend.hpp"
#include "core/precision.hpp"
#include "physics/mpas_atmphys_interface.hpp"
#include "physics/rrtmg/rrtmg_lw_tables.hpp"
#include "physics/rrtmg/rrtmg_sw_tables.hpp"

namespace cppmpas::phys {

// The DATA statements of the wrappers and of inirad / reicalc, device-copyable (they live with the RRTMG tables,
// physics/mpas_atmphys_radiation.hpp RrtmgTables): retab (LW 11547-11566), the standard profile pprof / tprof
// (LW 11726-11749), o3data's o3sum / ppsum / o3win / ppwin (LW 12345-12384).
struct RrtmgMpasConsts {
  real retab[95];
  real pprof[60], tprof[60];
  real o3sum[31], ppsum[31], o3win[31], ppwin[31];
};
void rrtmg_mpas_fill_consts(RrtmgMpasConsts& k);

// The capacity of the ozone levels of vinterp_ozn's column work (nOznLevels = 59 for MPAS's OZONE_PLEV.TBL)
inline constexpr int kMaxOznLevels = 64;

// The arguments of rrtmg_lwrad / rrtmg_swrad as MPAS's drivers pass them (mpas_atmphys_driver_radiation_lw.F:858-878,
// _sw.F:922-946): the mpas_atmphys_vars module arrays, (i) and (i, k) with ims = its = 1, kms = kts = 1, as kernel
// views (managed memory in the CUDA build).
struct RrtmgRadIn {
  int kte = 0;
  PCView p3d, p8w, pi3d, t3d, t8w, dz8w, qv3d, qc3d, qi3d, qs3d, cldfra3d;
  PCView o3clim;           // (i, 1:noznlevels)
  PCView pin;              // (1:noznlevels)
  int noznlevels = 0;
  PCView tsk, xland, xice, snow, xlat;  // (i)
  int julday = 0, icloud = 1, cldovrlp = 2, idcor = 0, o3input = 0;
  int has_reqc = 0, has_reqi = 0, has_reqs = 0;
};

struct RrtmgLwradOut {
  PView glw, olr, lwcf, lwupt, lwuptc, lwdnt, lwdntc, lwupb, lwupbc, lwdnb, lwdnbc;  // (i)
  PView rthratenlw, rre_cloud, rre_ice, rre_snow;                                    // (i, k)
};

// LW 11574-12305 on columns its..ite (its = 1). tables, k: device-accessible (RrtmgTables).
void rrtmg_lwrad(const RrtmgRadIn& in, PCView emiss, const RrtmgLwradOut& out, const rrtmg::RrtmgLwTables* tables,
                 const RrtmgMpasConsts* k, int its, int ite);

struct RrtmgSwradOut {
  PView coszr, gsw, swcf, swupt, swuptc, swdnt, swdntc, swupb, swupbc, swdnb, swdnbc, swddir, swddni, swddif;  // (i)
  PView rthratensw;                                                                                           // (i, k)
};

// SW 10025-10667 on columns its..ite (its = 1). aer_opt = 0 with MPAS's tauaer3d = 0, ssaaer3d = 1, asyaer3d = 0
// (radiation_sw_from_MPAS, mpas_atmphys_driver_radiation_sw.F:480-488): the values CPPWRF's spcvmc_sw takes.
void rrtmg_swrad(const RrtmgRadIn& in, PCView albedo, PCView xlong, real xtime, real gmt, real radt, real degrad,
                 real declin, real solcon, const RrtmgSwradOut& out, const rrtmg::RrtmgSwTables* tables,
                 const RrtmgMpasConsts* k, int its, int ite);

// vinterp_ozn (VI:16-176) for one column (ncol = pcols = 1, lchnk unused): o3vmr_out(1:pver) from pmid_in(1:pver),
// pin_in(1:levsiz), ozmix_in(1:levsiz); element k of each at P(k) / PIN(k) / O(k) / V(k) (1-based accessors).
// MPAS's copies into pin, ozmix and the flipped pmid (VI:73-86) and back (VI:160-165) are read in place: the same
// values. o3vmr_out(k) returns a reference. Its multi-column form (the LW driver's o3vmr diagnostic, mpas_atmphys_driver_radiation_lw.F:494-500) is
// column by column the same when every column's pmid is ordered (physics/mpas_atmphys_radiation.cpp).
template <class PM, class PI, class OZ, class OUT>
CPPMPAS_HD void vinterp_ozn_column(int pver, PM pmid_in, PI pin_in, int levsiz, OZ ozmix_in, OUT o3vmr_out) {
  auto PIN = [&](int k) { return pin_in(k); };
  auto O = [&](int k) { return ozmix_in(k); };
  auto P = [&](int k) { return pmid_in(pver - k + 1); };   // pmid(i, kk) = pmid_in(i, k), kk = pver - k + 1
  auto V = [&](int k, real x) { o3vmr_out(pver - k + 1) = x; };  // o3vmr_out(i, kk) = o3vmr(i, k)
  int kupper = 1;                    // VI:91-93
  for (int k = 1; k <= pver; ++k) {  // VI:95
    const int kkstart = levsiz < kupper ? levsiz : kupper;  // VI:100-103 (min over the one column)
    int kount = 0;
    bool found = false;
    for (int kk = kkstart; kk <= levsiz - 1; ++kk) {  // VI:109-131
      if (PIN(kk) < P(k) && P(k) <= PIN(kk + 1)) {
        kupper = kk;
        kount = kount + 1;
      }
      if (kount == 1) {  // VI:120-130 (kount == ncol)
        const real dpu = P(k) - PIN(kupper);
        const real dpl = PIN(kupper + 1) - P(k);
        V(k, (O(kupper) * dpl + O(kupper + 1) * dpu) / (dpl + dpu));
        found = true;
        break;  // goto 35
      }
    }
    if (found) continue;
    if (P(k) < PIN(1)) {  // VI:137-148
      V(k, O(1) * P(k) / PIN(1));
    } else if (P(k) > PIN(levsiz)) {
      V(k, O(levsiz));
    } else {
      const real dpu = P(k) - PIN(kupper);
      const real dpl = PIN(kupper + 1) - P(k);
      V(k, (O(kupper) * dpl + O(kupper + 1) * dpu) / (dpl + dpu));
    }
  }
}

}  // namespace cppmpas::phys
