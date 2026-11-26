// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// lsm (physics_wrf/module_sf_noahdrv.F:42-1771) and seaice_noah (physics_wrf/module_sf_noah_seaice_drv.F:18-517)
// of MPAS v8.4.2 with the arguments MPAS passes (mpas_atmphys_driver_lsm.F:705-749,
// mpas_atmphys_driver_seaice.F:488-516). The point work is CPPWRF's (physics/module_sf_noahdrv_kernels.hpp,
// physics/sf_noah.hpp, physics/sf_noah_ice.hpp); these functions run it for i = its..ite in lsm's order.
// Arrays are the driver's local (ims:ime, jms:jme) arrays (jms = jme = 1), soil arrays (ims:ime, 1:nsoil), the
// atmospheric ones (ims:ime, kms:kme); every pointer is the array's first element.
// Refused (hard errors, no verified port): ua_phys, fasdas /= 0, sf_urban_physics /= 0, rdlai2d, myj, frpcpn,
// opt_thcnd /= 1, num_soil_layers > noah::kMaxSoilLevels, kts /= kms; for seaice_noah seaice_albedo_opt,
// seaice_thickness_opt or seaice_snowdepth_opt /= 0.
#include "core/precision.hpp"
#include "physics/sf_noah.hpp"

namespace cppmpas::phys::noahdrv {

struct LsmArgs {
  const real *dz8w = nullptr, *p8w3d = nullptr, *t3d = nullptr, *qv3d = nullptr;  // (ims:ime, kms:kme)
  real* xland = nullptr;
  const real* xice = nullptr;
  const int *ivgtyp = nullptr, *isltyp = nullptr;
  const real *tmn = nullptr, *vegfra = nullptr, *shdmin = nullptr, *shdmax = nullptr, *snoalb = nullptr,
             *glw = nullptr, *swdown = nullptr, *rainbl = nullptr, *embck = nullptr, *sr = nullptr,
             *qgh = nullptr, *qz0 = nullptr, *chs = nullptr, *rib = nullptr, *cqs2 = nullptr;
  real *tsk = nullptr, *hfx = nullptr, *qfx = nullptr, *lh = nullptr, *grdflx = nullptr, *qsfc = nullptr,
       *chs2 = nullptr, *snow = nullptr, *snowc = nullptr, *snowh = nullptr, *canwat = nullptr, *smstav = nullptr,
       *smstot = nullptr, *sfcrunoff = nullptr, *udrunoff = nullptr, *acsnom = nullptr, *acsnow = nullptr,
       *snotime = nullptr, *snopcx = nullptr, *emiss = nullptr, *potevp = nullptr, *albedo = nullptr,
       *albbck = nullptr, *z0 = nullptr, *znt = nullptr, *lai = nullptr, *noahres = nullptr, *chklowq = nullptr;
  real *sh2o = nullptr, *smois = nullptr, *tslb = nullptr, *smcrel = nullptr;  // (ims:ime, 1:nsoil)
  const real* dzs = nullptr;                                                    // (1:nsoil)
  int isurban = 0, isice = 0;
  real dt = 0;
  bool myj = false, frpcpn = false, rdlai2d = false, ua_phys = false, usemonalb = true;
  int itimestep = 0, opt_thcnd = 1, fasdas = 0, num_soil_layers = 0, sf_urban_physics = 0;
  real xice_threshold = 0;
  int ims = 1, ime = 0, its = 1, ite = 0, kms = 1, kts = 1;
};
void lsm(const noah::NoahTables& tables, const LsmArgs& a);

struct SeaiceArgs {
  const real *dz8w = nullptr, *p8w3d = nullptr, *t3d = nullptr, *qv3d = nullptr;  // (ims:ime, kms:kme)
  const real *xice = nullptr, *snoalb2d = nullptr, *glw = nullptr, *swdown = nullptr, *rainbl = nullptr,
             *sr = nullptr, *qgh = nullptr, *rib = nullptr;
  real *tsk = nullptr, *hfx = nullptr, *qfx = nullptr, *lh = nullptr, *grdflx = nullptr, *qsfc = nullptr,
       *emiss = nullptr, *albedo = nullptr, *cqs2 = nullptr, *chs = nullptr, *chs2 = nullptr, *z02d = nullptr,
       *znt = nullptr, *snow = nullptr, *snowc = nullptr, *snowh2d = nullptr, *acsnow = nullptr, *acsnom = nullptr,
       *sfcrunoff = nullptr, *noahres = nullptr, *potevp = nullptr, *snopcx = nullptr;
  real* tslb = nullptr;  // (ims:ime, 1:nsoil)
  real dt = 0;
  bool frpcpn = false;
  int seaice_albedo_opt = 0, seaice_thickness_opt = 0, seaice_snowdepth_opt = 0;
  real seaice_albedo_default = 0, seaice_thickness_default = 0, seaice_snowdepth_max = 0, seaice_snowdepth_min = 0;
  real xice_threshold = 0;
  int num_soil_layers = 0, sf_urban_physics = 0;
  int ims = 1, ime = 0, its = 1, ite = 0, kms = 1, kts = 1;
};
void seaice_noah(const SeaiceArgs& a);

}  // namespace cppmpas::phys::noahdrv
