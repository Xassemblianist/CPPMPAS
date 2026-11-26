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
// Reference object code (bin/atmosphere_model, sha256 524a1cf4...): no FMA and no ymm register in any radiation
// routine. rrtmg_lwrad's only libm call is lround (NINT, LW 11822); rrtmg_swrad calls sincos twice and cos once
// (SW 10205-10206, the zenith angle; xxlat and declin take sincos) . In the kernels every libm call is scalar (exp,
// log, pow) except one libmvec lane: generate_stochastic_clouds' exponential-overlap alpha (_ZGVbN2v_exp, LW 2552 /
// SW 1860), reached only with config_radt_cld_overlap = 'exponential', which is refused here.
#include <vector>

#include "core/precision.hpp"
#include "physics/mpas_atmphys_interface.hpp"
#include "physics/rrtmg/rrtmg_lw_tables.hpp"
#include "physics/rrtmg/rrtmg_sw_tables.hpp"

namespace cppmpas::phys {

// module_ra_rrtmg_vinterp.F:16-176. pmid_in(pcols, pver), ozmix_in(pcols, levsiz), o3vmr_out(pcols, pver): Fortran
// order, column fastest; pin_in(levsiz).
void vinterp_ozn(int lchnk, int ncol, int pcols, int pver, const real* pmid_in, const real* pin_in, int levsiz,
                 const real* ozmix_in, real* o3vmr_out);

// inirad (LW 12308-12331) with o3data (LW 12334-12444): o3prof(kts:kte+1) from plev(kts:kte+2) (1-based arrays
// passed as pointers to their first element)
void inirad(real* o3prof, const real* plev, int kts, int kte);

// The arguments of rrtmg_lwrad / rrtmg_swrad as MPAS's drivers pass them (mpas_atmphys_driver_radiation_lw.F:858-878,
// _sw.F:922-946): the mpas_atmphys_vars module arrays, (i) and (i, k) with ims = its = 1, kms = kts = 1.
struct RrtmgRadIn {
  int kte = 0;
  const PArray *p3d = nullptr, *p8w = nullptr, *pi3d = nullptr, *t3d = nullptr, *t8w = nullptr, *dz8w = nullptr;
  const PArray *qv3d = nullptr, *qc3d = nullptr, *qi3d = nullptr, *qs3d = nullptr, *cldfra3d = nullptr;
  const PArray *re_cloud = nullptr, *re_ice = nullptr, *re_snow = nullptr;
  const PArray* o3clim = nullptr;  // (i, 1:noznlevels)
  const real* pin = nullptr;       // (1:noznlevels)
  int noznlevels = 0;
  const std::vector<real>*tsk = nullptr, *xland = nullptr, *xice = nullptr, *snow = nullptr, *xlat = nullptr;
  int julday = 0, icloud = 1, cldovrlp = 2, idcor = 0, o3input = 0;
  int has_reqc = 0, has_reqi = 0, has_reqs = 0;
};

struct RrtmgLwradOut {
  std::vector<real>*glw, *olr, *lwcf, *lwupt, *lwuptc, *lwdnt, *lwdntc, *lwupb, *lwupbc, *lwdnb, *lwdnbc;
  PArray *rthratenlw, *rre_cloud, *rre_ice, *rre_snow;
};

// LW 11574-12305 on columns its..ite (its = 1)
void rrtmg_lwrad(const RrtmgRadIn& in, const std::vector<real>& emiss, RrtmgLwradOut& out,
                 const rrtmg::RrtmgLwTables& tables, int its, int ite);

struct RrtmgSwradOut {
  std::vector<real>*coszr, *gsw, *swcf, *swupt, *swuptc, *swdnt, *swdntc, *swupb, *swupbc, *swdnb, *swdnbc;
  std::vector<real>*swddir, *swddni, *swddif;
  PArray* rthratensw;
};

// SW 10025-10667 on columns its..ite (its = 1). aer_opt = 0 with MPAS's tauaer3d = 0, ssaaer3d = 1, asyaer3d = 0
// (radiation_sw_from_MPAS, mpas_atmphys_driver_radiation_sw.F:480-488): the values CPPWRF's spcvmc_sw takes.
void rrtmg_swrad(const RrtmgRadIn& in, const std::vector<real>& albedo, const std::vector<real>& xlong, real xtime,
                 real gmt, real radt, real degrad, real declin, real solcon, RrtmgSwradOut& out,
                 const rrtmg::RrtmgSwTables& tables, int its, int ite);

}  // namespace cppmpas::phys
