// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// The point kernels of the Noah drivers: lsm of physics_wrf/module_sf_noahdrv.F (MPAS v8.4.2; KLsmPre,
// KLsmGlacial, KLsmPost) and seaice_noah of physics_wrf/module_sf_noah_seaice_drv.F (KSeaiceNoah), from CPPWRF
// (commit 7fac5bf, src/phys/driver/surface_driver.cpp lines 171-462 and 565-657, last changed in 8c2edeb). Shared
// with CPPWRF (SHARED_WITH_CPPWRF.md).
//
// Same code: MPAS's module_sf_noahdrv.F differs from WRF 4.8.0's phys/module_sf_noahdrv.F only in urban-canopy
// code (sf_urban_physics 1-3: the Chenghao Wang urban vegetation arrays, slucm_distributed_drag, lf_urb2d_s /
// z0_urb2d, FRC_URB2D INOUT for the mosaic LCZ fractions) and in lsm_mosaic (diff: 246 lines, none in lsm's
// sf_urban_physics = 0 path). MPAS's module_sf_noah_seaice_drv.F declares POTEVP, SNOPCX and NOAHRES OPTIONAL and
// updates them under PRESENT(); MPAS's driver passes all three (mpas_atmphys_driver_seaice.F:501-502), and the
// POTEVP accumulation moved after the QSFC update: independent statements, the same values.
// lsm's CYCLE, the land-sea test and the ICE = 1 / 0 / -1 cases are the kernels' returns and ice() codes; lsm runs
// them for one point after the other (module_sf_noahdrv.cpp). Loop 50 (ITIMESTEP = 1, module_sf_noahdrv.F:746-785)
// is per point, so it is KLsmPre's first block.
//
// Copied unchanged except CPPWRF_HD spelled CPPMPAS_HD, the namespace, and the adapters below:
//  * wrfc:: is module_sf_noahlsm's USE of mpas_atmphys_constants under -Dmpas (module_sf_noahlsm.F:3-4; WRF:
//    module_model_constants, the same values);
//  * rpow, libm::, raise_fatal as in physics/sf_noah_ice.hpp; kFatalNoahUrbanLai is CPPWRF's stop for the rdlai2d
//    path, which MPAS never takes (rdlai2d = .false., mpas_atmphys_driver_lsm.F:34).
// Object code of the reference executable: __module_sf_noahdrv_MOD_lsm calls pow@plt (APES, APELM) and exp@plt
// (E2SAT); __module_sf_noah_seaice_drv_MOD_seaice_noah the same; no libmvec call, no FMA, no ymm register.
#include "core/backend.hpp"
#include "core/libm.hpp"
#include "core/precision.hpp"
#include "physics/mpas_atmphys_constants.hpp"
#include "physics/sf_noah.hpp"
#include "physics/sf_noah_ice.hpp"

namespace cppmpas::phys::noahdrv {

using V2 = noah::View2<real>;
using C2 = noah::View2<const real>;
using V3 = noah::View3<real>;
using C3 = noah::View3<const real>;
using I2 = noah::View2<int>;
using CI2 = noah::View2<const int>;
constexpr int kMaxSoil = noah::kMaxSoilLevels;
namespace noah_ice = ::cppmpas::phys::noah_ice;
namespace noah = ::cppmpas::phys::noah;
using noah::raise_fatal;
constexpr int kFatalNoahUrbanLai = 16;
using noah_ice::rpow;
// module_sf_noahlsm.F:3-4 under -Dmpas: CP => cp, R_D => R_d, XLF => xlf, XLV => xlv, RHOWATER => rho_w,
// STBOLT => stbolt, KARMAN => karman of mpas_atmphys_constants
namespace wrfc {
constexpr real r_d = ::cppmpas::physics_constants::R_d;
constexpr real cp = ::cppmpas::physics_constants::cp;
constexpr real xlf = ::cppmpas::physics_constants::xlf;
constexpr real xlv = ::cppmpas::physics_constants::xlv;
constexpr real rhowater = ::cppmpas::physics_constants::rho_w;
constexpr real stbolt = ::cppmpas::physics_constants::stbolt;
}  // namespace wrfc

// ---- copied from CPPWRF 7fac5bf src/phys/driver/surface_driver.cpp:171-462 ----
// ---------------------------------------------------------------- Noah driver (lsm), before SFLX
struct LsmConst {
  // module_sf_noahdrv.F lsm parameters
  static constexpr real capa = wrfc::r_d / wrfc::cp;
  static constexpr real tresh = real(.95e0), a2 = real(17.67), a3 = real(273.15), a4 = real(29.65);
  static constexpr real a23m4 = a2 * (a3 - a4);
  static constexpr real row = real(1.e3), eliw = wrfc::xlf, rowliw = row * eliw;
};

struct KLsmPre {
  // inputs
  C3 dz8w, qv3d, p8w3d, t3d;
  C2 qgh, glw, swdown, rainbl, vegfra, shdmin, shdmax, snoalb, tmn, xice, chs, rib, sr, qz0;
  CI2 ivgtyp, isltyp;
  // in/out
  V2 tsk, albedo, albbck, z0, emiss, embck, snow, snowh, snowc, canwat, snotime, lai, smstav, smstot, chklowq,
      xland;
  V3 smois, tslb, sh2o, smcrel;
  // SFLX / SFLX_GLACIAL inputs (scratch)
  V2 lwdn, soldn, solnet, sfcprs, prcp, sfctmp, th2, q2k, q2sat, dqsdt2, ffrozp, shdfac, shmin, shmax, albbrd,
      snoalb1, tbot, z0brd, emissi, embrd, cmc, t1, snowhk, sneqv, chk, snotime1, ribb, xlai, soilw, sncovr;
  V3 stc, smc, swc;
  I2 vegtyp, soiltyp, slopetyp, ice, active;
  int nsoil, itimestep, isurban, isice, sf_urban_physics;
  int lcz[11];
  bool myj, frpcpn, rdlai2d;
  real dt, xice_threshold;
  // sf_urban_physics = 1: FRC_URB2D, TS_URB2D; the NATURAL category of VEGPARM.TBL and its SHDTBL
  C2 frc_urb, ts_urb;
  int natural;
  real shdtbl_natural;
  CPPMPAS_HD void operator()(int i, int j) const {
    using C = LsmConst;
    ice(i, j) = 0;
    active(i, j) = 0;
    if (itimestep == 1) {
      if ((xland(i, j) - real(1.5)) >= real(0.)) {
        smstav(i, j) = real(1.0);
        smstot(i, j) = real(1.0);
        for (int ns = 0; ns < nsoil; ++ns) {
          smois(i, j, ns) = real(1.0);
          tslb(i, j, ns) = real(273.16);
          smcrel(i, j, ns) = real(1.0);
        }
      } else if (xice(i, j) >= xice_threshold) {
        smstav(i, j) = real(1.0);
        smstot(i, j) = real(1.0);
        for (int ns = 0; ns < nsoil; ++ns) {
          smois(i, j, ns) = real(1.0);
          smcrel(i, j, ns) = real(1.0);
        }
      }
    }
    real psfc = p8w3d(i, j, 0);
    sfcprs(i, j) = (p8w3d(i, j, 1) + p8w3d(i, j, 0)) * real(0.5);
    q2k(i, j) = qv3d(i, j, 0) / (real(1.0) + qv3d(i, j, 0));
    soilw(i, j) = real(0.);
    real q2s = qgh(i, j) / (real(1.0) + qgh(i, j));
    if (myj && (q2k(i, j) >= q2s * C::tresh) && (q2k(i, j) < qz0(i, j))) {  // QZ0 test of MYJ (not used with sfclayrev)
      chklowq(i, j) = real(0.);
    } else {
      chklowq(i, j) = real(1.);
    }
    real sfct = t3d(i, j, 0);
    sfctmp(i, j) = sfct;
    real apes = rpow(real(1.e5) / psfc, C::capa);
    real apelm = rpow(real(1.e5) / sfcprs(i, j), C::capa);
    real sfcth2 = sfct * apelm;
    th2(i, j) = sfcth2 / apes;
    real emis = emiss(i, j);
    emissi(i, j) = emis;
    lwdn(i, j) = glw(i, j) * emis;
    real sold = swdown(i, j);
    soldn(i, j) = sold;
    solnet(i, j) = sold * (real(1.) - albedo(i, j));
    prcp(i, j) = rainbl(i, j) / dt;
    int vtyp = ivgtyp(i, j);
    int styp = isltyp(i, j);
    real shdf = vegfra(i, j) / real(100.);
    t1(i, j) = tsk(i, j);
    chk(i, j) = chs(i, j);
    shmin(i, j) = shdmin(i, j) / real(100.);
    shmax(i, j) = shdmax(i, j) / real(100.);
    real sneq = snow(i, j) * real(0.001);
    real snwh = snowh(i, j);
    sncovr(i, j) = snowc(i, j);
    if (frpcpn) {
      ffrozp(i, j) = sr(i, j);
    } else {
      ffrozp(i, j) = sfct <= real(273.15) ? real(1.0) : real(0.0);
    }
    slopetyp(i, j) = 1;
    if ((xland(i, j) - real(1.5)) >= real(0.)) {  // open water: no land-surface model
      vegtyp(i, j) = vtyp;
      soiltyp(i, j) = styp;
      return;
    }
    int ic;
    if (xice(i, j) >= xice_threshold) {
      ic = 1;
    } else if (vtyp == isice) {
      ic = -1;
    } else {
      ic = 0;
    }
    real dq = q2s * C::a23m4 / ((sfct - C::a4) * (sfct - C::a4));
    if (snow(i, j) > real(0.0)) {
      real sfctsno = sfct;
      real e2sat = real(611.2) * libm::exp(real(6174.) * (real(1.) / real(273.15) - real(1.) / sfctsno));
      real q2sati = real(0.622) * e2sat / (sfcprs(i, j) - e2sat);
      q2sati = q2sati / (real(1.0) + q2sati);
      if (t1(i, j) > real(273.14)) {
        q2s = q2s * (real(1.) - snowc(i, j)) + q2sati * snowc(i, j);
        dq = dq * (real(1.) - snowc(i, j)) + q2sati * real(6174.) / (sfctsno * sfctsno) * snowc(i, j);
      } else {
        q2s = q2sati;
        dq = q2sati * real(6174.) / (sfctsno * sfctsno);
      }
      if (t1(i, j) > real(273.) && snowc(i, j) > real(0.) && swdown(i, j) > real(10.)) dq = dq * (real(1.) - snowc(i, j));
    }
    q2sat(i, j) = q2s;
    dqsdt2(i, j) = dq;
    tbot(i, j) = tmn(i, j);
    if (isurban == 1) {  // USGS land data
      if (vtyp == 25) shdf = real(0.0000);
      if (vtyp == 26) shdf = real(0.0000);
      if (vtyp == 27) shdf = real(0.0000);
    }
    if (styp == 14 && xice(i, j) == real(0.)) styp = 7;
    snoalb1(i, j) = snoalb(i, j);
    cmc(i, j) = canwat(i, j) / real(1000.);
    albbrd(i, j) = albbck(i, j);
    z0brd(i, j) = z0(i, j);
    embrd(i, j) = embck(i, j);
    snotime1(i, j) = snotime(i, j);
    ribb(i, j) = rib(i, j);
    for (int ns = 0; ns < nsoil; ++ns) {
      smc(i, j, ns) = smois(i, j, ns);
      stc(i, j, ns) = tslb(i, j, ns);
      swc(i, j, ns) = sh2o(i, j, ns);
    }
    if ((sneq != real(0.) && snwh == real(0.)) || (snwh <= sneq)) snwh = real(5.) * sneq;
    // sf_urban_physics = 0: urban (and LCZ) points use the urban vegetation category; 1 (UCM): the natural surface
    // of the urban point is VEGPARM.TBL's NATURAL category, SFLX sees the rural skin temperature
    bool urb = vtyp == isurban;
    for (int n = 0; n < 11; ++n) urb = urb || vtyp == lcz[n];
    if (sf_urban_physics == 0 && urb) vtyp = isurban;
    if (sf_urban_physics == 1 && urb) {
      vtyp = natural;
      shdf = shdtbl_natural;
      albbrd(i, j) = real(0.2);
      emissi(i, j) = real(0.98);
      lwdn(i, j) = glw(i, j) * emissi(i, j);
      solnet(i, j) = sold * (real(1.0) - real(0.2));  // SOLDN * (1.0 - ALBEDOK), ALBEDOK = 0.2
      if (frc_urb(i, j) < real(0.99)) {
        t1(i, j) = (tsk(i, j) - frc_urb(i, j) * ts_urb(i, j)) / (real(1) - frc_urb(i, j));
      } else {
        t1(i, j) = tsk(i, j);
      }
    }
    if (rdlai2d) {
      if (shdf > real(0.0) && lai(i, j) <= real(0.0)) lai(i, j) = real(0.01);
      xlai(i, j) = lai(i, j);
      // SFLX turns SHDFAC 0 into 0.05 here, skips CANRES (XLAI = 0) and TRANSP reads the stale PC (R17 U1)
      if (ic == 0 && vtyp == isurban && shdf <= real(0.0) && lai(i, j) <= real(0.0)) raise_fatal(kFatalNoahUrbanLai);
    }
    vegtyp(i, j) = vtyp;
    soiltyp(i, j) = styp;
    shdfac(i, j) = shdf;
    sneqv(i, j) = sneq;
    snowhk(i, j) = snwh;
    ice(i, j) = ic;
    if (ic == 1) {  // sea ice: handled by seaice_noah
      for (int ns = 0; ns < nsoil; ++ns) sh2o(i, j, ns) = real(1.0);
      lai(i, j) = real(0.01);
      ice(i, j) = 1;
      return;
    }
    if (ic == 0) active(i, j) = 1;
  }
};

// SFLX_GLACIAL for land-ice points (ICE == -1); outputs in the SFLX scratch of the driver.
struct KLsmGlacial {
  CI2 ice;
  V2 lwdn, solnet, sfcprs, prcp, sfctmp, th2, q2k, q2sat, dqsdt2, ffrozp, albbrd, snoalb1, tbot, z0brd, emissi,
      embrd, t1, snowhk, sneqv, chk, snotime1, ribb, xlai, soilw, soilm, runoff1, runoff2, runoff3, eta, sheat,
      eta_kin, etp, ssoil, flx1, flx2, flx3, snomlt, sncovr, q1, albedok, z0k;
  V3 stc, smc, swc, smav;
  real dt, lvcoef;
  real sldpth[kMaxSoil];
  int nsoil;
  CPPMPAS_HD void operator()(int i, int j) const {
    if (ice(i, j) != -1) return;
    soilm(i, j) = real(0.0);
    xlai(i, j) = real(0.01);
    runoff2(i, j) = real(0.0);
    runoff3(i, j) = real(0.0);
    for (int ns = 0; ns < nsoil; ++ns) {
      swc(i, j, ns) = real(1.0);
      smc(i, j, ns) = real(1.0);
      smav(i, j, ns) = real(1.0);
    }
    real stcl[kMaxSoil];
    for (int ns = 0; ns < nsoil; ++ns) stcl[ns] = stc(i, j, ns);
    real em = emissi(i, j), tt = t1(i, j), sh = snowhk(i, j), se = sneqv(i, j), sc = sncovr(i, j),
         st = snotime1(i, j), rb = ribb(i, j);
    noah_ice::glacial::Out o{};
    noah_ice::glacial::sflx_glacial(ffrozp(i, j), dt, nsoil, sldpth, lwdn(i, j), solnet(i, j), sfcprs(i, j),
                                    prcp(i, j), sfctmp(i, j), q2k(i, j), th2(i, j), q2sat(i, j), dqsdt2(i, j),
                                    albbrd(i, j), snoalb1(i, j), tbot(i, j), z0brd(i, j), em, embrd(i, j), tt, stcl,
                                    sh, se, chk(i, j), sc, st, rb, lvcoef, o);
    for (int ns = 0; ns < nsoil; ++ns) stc(i, j, ns) = stcl[ns];
    emissi(i, j) = em;
    t1(i, j) = tt;
    snowhk(i, j) = sh;
    sneqv(i, j) = se;
    sncovr(i, j) = sc;
    snotime1(i, j) = st;
    ribb(i, j) = rb;
    albedok(i, j) = o.albedo;
    z0k(i, j) = o.z0;
    eta(i, j) = o.eta;
    sheat(i, j) = o.sheat;
    eta_kin(i, j) = o.eta_kinematic;
    etp(i, j) = o.etp;
    ssoil(i, j) = o.ssoil;
    flx1(i, j) = o.flx1;
    flx2(i, j) = o.flx2;
    flx3(i, j) = o.flx3;
    snomlt(i, j) = o.snomlt;
    runoff1(i, j) = o.runoff1;
    q1(i, j) = o.q1;
  }
};

// lsm after SFLX / SFLX_GLACIAL (land points, ICE /= 1)
struct KLsmPost {
  CI2 ice;
  C2 xland;
  C2 lwdn, solnet, prcp, ffrozp, cmc, sneqv, snowhk, albedok, albbrd, z0brd, emissi, z0k, t1, sheat, etp, eta_kin,
      eta, ssoil, sncovr, snotime1, q1, flx1, flx2, flx3, soilw, soilm, runoff1, runoff2, snomlt, xlai, cqs2;
  C3 smc, stc, swc, smav;
  V2 lai, canwat, snow, snowh, albedo, albbck, z0, emiss, znt, tsk, hfx, potevp, qfx, lh, grdflx, snowc, chs2,
      snotime, qsfc, noahres, smstav, smstot, sfcrunoff, udrunoff, acsnow, acsnom, snopcx;
  V3 smois, tslb, sh2o, smcrel;
  real dt;
  int nsoil;
  CPPMPAS_HD void operator()(int i, int j) const {
    using C = LsmConst;
    if ((xland(i, j) - real(1.5)) >= real(0.) || ice(i, j) == 1) return;
    const real fdtliw = dt / C::rowliw;
    const real fdtw = dt / (wrfc::xlv * wrfc::rhowater);
    lai(i, j) = xlai(i, j);
    canwat(i, j) = cmc(i, j) * real(1000.);
    snow(i, j) = sneqv(i, j) * real(1000.);
    snowh(i, j) = snowhk(i, j);
    albedo(i, j) = albedok(i, j);
    albbck(i, j) = albbrd(i, j);
    z0(i, j) = z0brd(i, j);
    emiss(i, j) = emissi(i, j);
    znt(i, j) = z0k(i, j);
    real tt = t1(i, j);
    tsk(i, j) = tt;
    hfx(i, j) = sheat(i, j);
    potevp(i, j) = potevp(i, j) + etp(i, j) * fdtw;
    qfx(i, j) = eta_kin(i, j);
    lh(i, j) = eta(i, j);
    grdflx(i, j) = ssoil(i, j);
    snowc(i, j) = sncovr(i, j);
    chs2(i, j) = cqs2(i, j);
    snotime(i, j) = snotime1(i, j);
    qsfc(i, j) = q1(i, j) / (real(1.0) - q1(i, j));
    for (int ns = 0; ns < nsoil; ++ns) {
      smois(i, j, ns) = smc(i, j, ns);
      tslb(i, j, ns) = stc(i, j, ns);
      sh2o(i, j, ns) = swc(i, j, ns);
    }
    noahres(i, j) = (solnet(i, j) + lwdn(i, j)) - sheat(i, j) + ssoil(i, j) - eta(i, j) -
                    (emissi(i, j) * wrfc::stbolt * ((tt * tt) * (tt * tt))) - flx1(i, j) - flx2(i, j) - flx3(i, j);
    smstav(i, j) = soilw(i, j);
    smstot(i, j) = soilm(i, j) * real(1000.);
    for (int ns = 0; ns < nsoil; ++ns) smcrel(i, j, ns) = smav(i, j, ns);
    sfcrunoff(i, j) = sfcrunoff(i, j) + runoff1(i, j) * dt * real(1000.0);
    udrunoff(i, j) = udrunoff(i, j) + runoff2(i, j) * dt * real(1000.0);
    if (ffrozp(i, j) > real(0.5)) acsnow(i, j) = acsnow(i, j) + prcp(i, j) * dt;
    if (snow(i, j) > real(0.)) {
      acsnom(i, j) = acsnom(i, j) + snomlt(i, j) * real(1000.);
      snopcx(i, j) = snopcx(i, j) - snomlt(i, j) / fdtliw;
    }
  }
};

// ---- copied from CPPWRF 7fac5bf src/phys/driver/surface_driver.cpp:565-657 ----
// module_sf_noah_seaice_drv.F seaice_noah
struct KSeaiceNoah {
  C3 t3d, qv3d, p8w3d, dz8w;
  C2 glw, swdown, rainbl, snoalb2d, qgh, xice, rib, sr;
  V2 emiss, albedo, z02d, tsk, snow, snowc, snowh2d, chs, chs2, cqs2, znt, lh, hfx, qfx, potevp, grdflx, qsfc,
      acsnow, acsnom, snopcx, sfcrunoff, noahres;
  V3 tslb;
  noah_ice::seaice::Opts op;
  int thickness_opt;
  real thickness_default;
  real dt, xice_threshold;
  int nsoil;
  bool frpcpn;
  CPPMPAS_HD void operator()(int i, int j) const {
    using C = LsmConst;
    const real fdtliw = dt / C::rowliw;
    const real fdtw = dt / (wrfc::xlv * wrfc::rhowater);
    if (xice(i, j) < xice_threshold) return;  // (ICEDEPTH/SNOWSI resets only for opt 1)
    real sithick = thickness_default;  // SEAICE_THICKNESS_OPT = 0 (opt 1 needs ICEDEPTH input)
    real sfctmp = t3d(i, j, 0);
    real t1 = tsk(i, j);
    real snoalb = snoalb2d(i, j);
    real emissi = emiss(i, j);
    real lwdn = glw(i, j) * emissi;
    real sneqv = snow(i, j) * real(0.001);
    real snowh = snowh2d(i, j);
    real sncovr = snowc(i, j);
    real solnet = swdown(i, j) * (real(1.) - albedo(i, j));
    real sfcprs = (p8w3d(i, j, 1) + p8w3d(i, j, 0)) * real(0.5);
    real psfc = p8w3d(i, j, 0);
    real q2 = qv3d(i, j, 0) / (real(1.0) + qv3d(i, j, 0));
    real apes = rpow(real(1.e5) / psfc, C::capa);
    real apelm = rpow(real(1.e5) / sfcprs, C::capa);
    real th2 = (sfctmp * apelm) / apes;
    real q2sat = qgh(i, j) / (real(1.0) + qgh(i, j));
    real dqsdt2 = q2sat * C::a23m4 / ((sfctmp - C::a4) * (sfctmp - C::a4));
    if (snow(i, j) > real(0.0)) {
      real sfctsno = sfctmp;
      real e2sat = real(611.2) * libm::exp(real(6174.) * (real(1.) / real(273.15) - real(1.) / sfctsno));
      real q2sati = real(0.622) * e2sat / (sfcprs - e2sat);
      q2sati = q2sati / (real(1.0) + q2sati);
      if (t1 > real(273.14)) {
        q2sat = q2sat * (real(1.) - snowc(i, j)) + q2sati * snowc(i, j);
        dqsdt2 = dqsdt2 * (real(1.) - snowc(i, j)) + q2sati * real(6174.) / (sfctsno * sfctsno) * snowc(i, j);
      } else {
        q2sat = q2sati;
        dqsdt2 = q2sati * real(6174.) / (sfctsno * sfctsno);
      }
      if ((t1 > real(273.)) && (snowc(i, j) > real(0.0))) dqsdt2 = dqsdt2 * (real(1.) - snowc(i, j));
    }
    real prcp = rainbl(i, j) / dt;
    real ffrozp;
    if (frpcpn) {
      ffrozp = sr(i, j);
    } else {
      ffrozp = sfctmp <= real(273.15) ? real(1.0) : real(0.0);
    }
    real tbot = real(271.36);
    real z0brd = z02d(i, j);
    real stc[kMaxSoil];
    for (int ns = 0; ns < nsoil; ++ns) stc[ns] = tslb(i, j, ns);
    real ch = chs(i, j);
    real ribb = rib(i, j);
    noah_ice::seaice::Out o{};
    noah_ice::seaice::sflx_seaice(op, ffrozp, dt, nsoil, sithick, lwdn, solnet, sfcprs, prcp, sfctmp, q2, th2, q2sat,
                                  dqsdt2, snoalb, tbot, z0brd, emissi, t1, stc, snowh, sneqv, ch, sncovr,
                                  ribb, o);
    albedo(i, j) = o.albedo;
    emiss(i, j) = emissi;
    tsk(i, j) = t1;
    z02d(i, j) = z0brd;
    snowh2d(i, j) = snowh;
    snowc(i, j) = sncovr;
    snow(i, j) = sneqv * real(1000.);
    for (int ns = 0; ns < nsoil; ++ns) tslb(i, j, ns) = stc[ns];
    znt(i, j) = o.z0;
    lh(i, j) = o.eta;
    hfx(i, j) = o.sheat;
    qfx(i, j) = o.eta_kinematic;
    potevp(i, j) = potevp(i, j) + o.etp * fdtw;
    grdflx(i, j) = o.ssoil;
    chs2(i, j) = cqs2(i, j);
    if (o.q1 > qsfc(i, j)) cqs2(i, j) = chs(i, j);
    qsfc(i, j) = o.q1 / (real(1.0) - o.q1);
    if (ffrozp > real(0.5)) acsnow(i, j) = acsnow(i, j) + prcp * dt;
    acsnom(i, j) = acsnom(i, j) + o.snomlt * real(1000.);
    snopcx(i, j) = snopcx(i, j) - o.snomlt / fdtliw;
    sfcrunoff(i, j) = sfcrunoff(i, j) + o.runoff1 * dt * real(1000.0);
    noahres(i, j) = (solnet + lwdn) - o.sheat + o.ssoil - o.eta - (emissi * wrfc::stbolt * ((t1 * t1) * (t1 * t1))) -
                    o.flx1 - o.flx2 - o.flx3;
    (void)thickness_opt;
  }
};

}  // namespace cppmpas::phys::noahdrv
