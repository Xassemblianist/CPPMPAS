// SPDX-License-Identifier: GPL-3.0-or-later
// lsm (physics_wrf/module_sf_noahdrv.F:42-1771) and seaice_noah (module_sf_noah_seaice_drv.F:18-517) with MPAS's
// arguments; see module_sf_noahdrv.hpp.
#include "physics/module_sf_noahdrv.hpp"

#include <string>
#include <vector>

#include "core/error.hpp"
#include "physics/module_sf_noahdrv_kernels.hpp"
#include "physics/mpas_atmphys_interface.hpp"

namespace cppmpas::phys::noahdrv {

namespace {

// A driver array (ims:ime, jms:jme) / (ims:ime, 1:n, jms:jme) as CPPWRF's View2 / View3 with i 0-based and j = 0
template <class T> noah::View2<T> v2(T* p) { return noah::View2<T>{p, 0}; }
template <class T> noah::View3<T> v3(T* p, int ni) { return noah::View3<T>{p, 0, ni}; }

// The locals of lsm that live across SFLX / SFLX_GLACIAL (module_sf_noahdrv.F:415-445), one per point of
// ims..ime: CPPWRF's kernels pass them between KLsmPre, the SFLX kernel, KLsmGlacial and KLsmPost by point.
struct LsmLocals {
  int ni, ns;
  WorkArray lwdn, soldn, solnet, sfcprs, prcp, sfctmp, th2, q2k, q2sat, dqsdt2, ffrozp, shdfac, shmin, shmax,
      albbrd, snoalb1, tbot, z0brd, emissi, embrd, cmc, t1, snowhk, sneqv, chk, snotime1, ribb, xlai, soilw, sncovr,
      albedok, cm, z0k, eta, sheat, eta_kin, fdown, ec, edir, ett, esnow, drip, dew, beta, etp, ssoil, flx1, flx2,
      flx3, snomlt, runoff1, runoff2, runoff3, rc, pc, rsmin, soilm, q1;
  WorkArray stc, smc, swc, et, smav;
  WorkArrayInt vegtyp, soiltyp, slopetyp, ice, active;
  LsmLocals(int ni_, int ns_) : ni(ni_), ns(ns_) {
    for (auto* f : {&lwdn, &soldn, &solnet, &sfcprs, &prcp, &sfctmp, &th2, &q2k, &q2sat, &dqsdt2, &ffrozp, &shdfac,
                    &shmin, &shmax, &albbrd, &snoalb1, &tbot, &z0brd, &emissi, &embrd, &cmc, &t1, &snowhk, &sneqv,
                    &chk, &snotime1, &ribb, &xlai, &soilw, &sncovr, &albedok, &cm, &z0k, &eta, &sheat, &eta_kin,
                    &fdown, &ec, &edir, &ett, &esnow, &drip, &dew, &beta, &etp, &ssoil, &flx1, &flx2, &flx3, &snomlt,
                    &runoff1, &runoff2, &runoff3, &rc, &pc, &rsmin, &soilm, &q1})
      f->assign(ni, real(0));
    for (auto* f : {&stc, &smc, &swc, &et, &smav}) f->assign(std::size_t(ni) * ns, real(0));
    for (auto* f : {&vegtyp, &soiltyp, &slopetyp, &ice, &active}) f->assign(ni, 0);
  }
  noah::View2<real> w(WorkArray& a) { return v2(a.data()); }
  noah::View2<const real> r(WorkArray& a) { return v2<const real>(a.data()); }
  noah::View3<real> w3(WorkArray& a) { return v3(a.data(), ni); }
  noah::View3<const real> r3(WorkArray& a) { return v3<const real>(a.data(), ni); }
  noah::View2<int> wi(WorkArrayInt& a) { return v2(a.data()); }
  noah::View2<const int> ri(WorkArrayInt& a) { return v2<const int>(a.data()); }
};

void refuse(bool bad, const char* what) {
  if (bad) throw Error(std::string("module_sf_noahdrv.F lsm: ") + what + " is not ported");
}

}  // namespace

// module_sf_noahdrv.F:42-1771, sf_urban_physics = 0, for the points its..ite (j = jts = jte)
void lsm(const noah::NoahTables& tables, const LsmArgs& a) {
  refuse(a.ua_phys, "ua_phys = .true. (UA Noah snow physics)");
  refuse(a.fasdas != 0, "fasdas /= 0");
  refuse(a.sf_urban_physics != 0, "sf_urban_physics /= 0 (urban canopy models)");
  refuse(a.rdlai2d, "rdlai2d = .true.");
  refuse(a.myj, "myj = .true.");
  refuse(a.frpcpn, "frpcpn = .true.");
  refuse(a.opt_thcnd != 1, "opt_thcnd /= 1");
  refuse(a.kts != a.kms, "kts /= kms");
  // one layer: SHFLX / SMFLX / the ice and glacial HRT read layer 2 of zsoil, stc, smc (undefined in MPAS)
  if (a.num_soil_layers < 2 || a.num_soil_layers > noah::kMaxSoilLevels)
    throw Error("module_sf_noahdrv.F lsm: num_soil_layers = " + std::to_string(a.num_soil_layers) +
                " is outside 2.." + std::to_string(noah::kMaxSoilLevels));
  const int ni = a.ime - a.ims + 1, nsoil = a.num_soil_layers;
  LsmLocals L(ni, nsoil);
  auto c2 = [](const real* p) { return v2<const real>(p); };
  auto w2 = [](real* p) { return v2<real>(p); };

  KLsmPre pre{};
  pre.dz8w = v3<const real>(a.dz8w, ni); pre.qv3d = v3<const real>(a.qv3d, ni);
  pre.p8w3d = v3<const real>(a.p8w3d, ni); pre.t3d = v3<const real>(a.t3d, ni);
  pre.qgh = c2(a.qgh); pre.glw = c2(a.glw); pre.swdown = c2(a.swdown); pre.rainbl = c2(a.rainbl);
  pre.vegfra = c2(a.vegfra); pre.shdmin = c2(a.shdmin); pre.shdmax = c2(a.shdmax); pre.snoalb = c2(a.snoalb);
  pre.tmn = c2(a.tmn); pre.xice = c2(a.xice); pre.chs = c2(a.chs); pre.rib = c2(a.rib); pre.sr = c2(a.sr);
  pre.qz0 = c2(a.qz0);
  pre.ivgtyp = v2<const int>(a.ivgtyp); pre.isltyp = v2<const int>(a.isltyp);
  pre.tsk = w2(a.tsk); pre.albedo = w2(a.albedo); pre.albbck = w2(a.albbck); pre.z0 = w2(a.z0);
  pre.emiss = w2(a.emiss); pre.embck = v2<real>(const_cast<real*>(a.embck)); pre.snow = w2(a.snow);
  pre.snowh = w2(a.snowh); pre.snowc = w2(a.snowc); pre.canwat = w2(a.canwat); pre.snotime = w2(a.snotime);
  pre.lai = w2(a.lai); pre.smstav = w2(a.smstav); pre.smstot = w2(a.smstot); pre.chklowq = w2(a.chklowq);
  pre.xland = w2(a.xland);
  pre.smois = v3(a.smois, ni); pre.tslb = v3(a.tslb, ni); pre.sh2o = v3(a.sh2o, ni); pre.smcrel = v3(a.smcrel, ni);
  pre.lwdn = L.w(L.lwdn); pre.soldn = L.w(L.soldn); pre.solnet = L.w(L.solnet); pre.sfcprs = L.w(L.sfcprs);
  pre.prcp = L.w(L.prcp); pre.sfctmp = L.w(L.sfctmp); pre.th2 = L.w(L.th2); pre.q2k = L.w(L.q2k);
  pre.q2sat = L.w(L.q2sat); pre.dqsdt2 = L.w(L.dqsdt2); pre.ffrozp = L.w(L.ffrozp); pre.shdfac = L.w(L.shdfac);
  pre.shmin = L.w(L.shmin); pre.shmax = L.w(L.shmax); pre.albbrd = L.w(L.albbrd); pre.snoalb1 = L.w(L.snoalb1);
  pre.tbot = L.w(L.tbot); pre.z0brd = L.w(L.z0brd); pre.emissi = L.w(L.emissi); pre.embrd = L.w(L.embrd);
  pre.cmc = L.w(L.cmc); pre.t1 = L.w(L.t1); pre.snowhk = L.w(L.snowhk); pre.sneqv = L.w(L.sneqv);
  pre.chk = L.w(L.chk); pre.snotime1 = L.w(L.snotime1); pre.ribb = L.w(L.ribb); pre.xlai = L.w(L.xlai);
  pre.soilw = L.w(L.soilw); pre.sncovr = L.w(L.sncovr);
  pre.stc = L.w3(L.stc); pre.smc = L.w3(L.smc); pre.swc = L.w3(L.swc);
  pre.vegtyp = L.wi(L.vegtyp); pre.soiltyp = L.wi(L.soiltyp); pre.slopetyp = L.wi(L.slopetyp);
  pre.ice = L.wi(L.ice); pre.active = L.wi(L.active);
  pre.nsoil = nsoil; pre.itimestep = a.itimestep; pre.isurban = a.isurban; pre.isice = a.isice;
  pre.sf_urban_physics = a.sf_urban_physics;
  for (int n = 0; n < 11; ++n) pre.lcz[n] = tables.lcz[n];
  pre.natural = 0; pre.shdtbl_natural = 0;  // sf_urban_physics = 0: not read
  pre.myj = a.myj; pre.frpcpn = a.frpcpn; pre.rdlai2d = a.rdlai2d; pre.dt = a.dt;
  pre.xice_threshold = a.xice_threshold;

  // SFLX for ICE = 0 (module_sf_noahdrv.F:1102-1129)
  noah::NoahKernel nk{};
  nk.lwdn = L.r(L.lwdn); nk.soldn = L.r(L.soldn); nk.solnet = L.r(L.solnet); nk.sfcprs = L.r(L.sfcprs);
  nk.prcp = L.r(L.prcp); nk.sfctmp = L.r(L.sfctmp); nk.th2 = L.r(L.th2); nk.q2 = L.r(L.q2k);
  nk.sfcspd = L.r(L.q2k);  // SFCSPD: DUMMY in lsm (not used by SFLX)
  nk.q2sat = L.r(L.q2sat); nk.dqsdt2 = L.r(L.dqsdt2); nk.ffrozp = L.r(L.ffrozp);
  nk.vegtyp = L.ri(L.vegtyp); nk.soiltyp = L.ri(L.soiltyp); nk.slopetyp = L.ri(L.slopetyp);
  nk.shdfac = L.r(L.shdfac); nk.shdmin = L.r(L.shmin); nk.shdmax = L.r(L.shmax);
  nk.alb = L.r(L.albbrd); nk.snoalb = L.r(L.snoalb1); nk.tbot = L.r(L.tbot);
  nk.z0brd = L.w(L.z0brd); nk.embrd = L.w(L.embrd);
  nk.cmc = L.w(L.cmc); nk.t1 = L.w(L.t1); nk.stc = L.w3(L.stc); nk.smc = L.w3(L.smc); nk.sh2o = L.w3(L.swc);
  nk.snowh = L.w(L.snowhk); nk.sneqv = L.w(L.sneqv); nk.albedo = L.w(L.albedok); nk.ch = L.w(L.chk);
  nk.cm = L.w(L.cm); nk.emissi = L.w(L.emissi); nk.z0 = L.w(L.z0k); nk.xlai = L.w(L.xlai);
  nk.snotime1 = L.w(L.snotime1); nk.ribb = L.w(L.ribb);
  nk.eta = L.w(L.eta); nk.sheat = L.w(L.sheat); nk.eta_kinematic = L.w(L.eta_kin); nk.fdown = L.w(L.fdown);
  nk.ec = L.w(L.ec); nk.edir = L.w(L.edir); nk.et = L.w3(L.et); nk.ett = L.w(L.ett); nk.esnow = L.w(L.esnow);
  nk.drip = L.w(L.drip); nk.dew = L.w(L.dew); nk.beta = L.w(L.beta); nk.etp = L.w(L.etp);
  nk.ssoil = L.w(L.ssoil); nk.flx1 = L.w(L.flx1); nk.flx2 = L.w(L.flx2); nk.flx3 = L.w(L.flx3);
  nk.snomlt = L.w(L.snomlt); nk.sncovr = L.w(L.sncovr); nk.runoff1 = L.w(L.runoff1);
  nk.runoff2 = L.w(L.runoff2); nk.runoff3 = L.w(L.runoff3); nk.rc = L.w(L.rc); nk.pc = L.w(L.pc);
  nk.rsmin = L.w(L.rsmin); nk.soilw = L.w(L.soilw); nk.soilm = L.w(L.soilm); nk.q1 = L.w(L.q1);
  nk.smav = L.w3(L.smav);
  nk.dt = a.dt; nk.zlvl = real(0);  // ZLVL = 0.5 DZ8W(i,1,j) (:815) is not read by SFLX
  nk.nsoil = nsoil;
  for (int k = 0; k < nsoil; ++k) nk.sldpth[k] = a.dzs[k];  // SLDPTH(NS) = DZS(NS) (:740-742)
  nk.isurban = a.isurban;
  KernelCopy<noah::NoahTables> tables_k;
  nk.tables = kernel_ptr(tables, tables_k);
  nk.rdlai2d = a.rdlai2d; nk.usemonalb = a.usemonalb; nk.opt_thcnd = a.opt_thcnd;
  nk.active = L.ri(L.active); nk.alb_out = L.w(L.albbrd);

  // SFLX_GLACIAL for ICE = -1 (:1137-1180)
  KLsmGlacial gl{};
  gl.ice = L.ri(L.ice); gl.lwdn = L.w(L.lwdn); gl.solnet = L.w(L.solnet); gl.sfcprs = L.w(L.sfcprs);
  gl.prcp = L.w(L.prcp); gl.sfctmp = L.w(L.sfctmp); gl.th2 = L.w(L.th2); gl.q2k = L.w(L.q2k);
  gl.q2sat = L.w(L.q2sat); gl.dqsdt2 = L.w(L.dqsdt2); gl.ffrozp = L.w(L.ffrozp); gl.albbrd = L.w(L.albbrd);
  gl.snoalb1 = L.w(L.snoalb1); gl.tbot = L.w(L.tbot); gl.z0brd = L.w(L.z0brd); gl.emissi = L.w(L.emissi);
  gl.embrd = L.w(L.embrd); gl.t1 = L.w(L.t1); gl.snowhk = L.w(L.snowhk); gl.sneqv = L.w(L.sneqv);
  gl.chk = L.w(L.chk); gl.snotime1 = L.w(L.snotime1); gl.ribb = L.w(L.ribb); gl.xlai = L.w(L.xlai);
  gl.soilw = L.w(L.soilw); gl.soilm = L.w(L.soilm); gl.runoff1 = L.w(L.runoff1); gl.runoff2 = L.w(L.runoff2);
  gl.runoff3 = L.w(L.runoff3); gl.eta = L.w(L.eta); gl.sheat = L.w(L.sheat); gl.eta_kin = L.w(L.eta_kin);
  gl.etp = L.w(L.etp); gl.ssoil = L.w(L.ssoil); gl.flx1 = L.w(L.flx1); gl.flx2 = L.w(L.flx2);
  gl.flx3 = L.w(L.flx3); gl.snomlt = L.w(L.snomlt); gl.sncovr = L.w(L.sncovr); gl.q1 = L.w(L.q1);
  gl.albedok = L.w(L.albedok); gl.z0k = L.w(L.z0k);
  gl.stc = L.w3(L.stc); gl.smc = L.w3(L.smc); gl.swc = L.w3(L.swc); gl.smav = L.w3(L.smav);
  gl.dt = a.dt; gl.lvcoef = tables.lvcoef_data; gl.nsoil = nsoil;
  for (int k = 0; k < nsoil; ++k) gl.sldpth[k] = a.dzs[k];

  // after SFLX / SFLX_GLACIAL (:1184-1311, 1569-1588)
  KLsmPost post{};
  post.ice = L.ri(L.ice); post.xland = c2(a.xland);
  post.lwdn = L.r(L.lwdn); post.solnet = L.r(L.solnet); post.prcp = L.r(L.prcp); post.ffrozp = L.r(L.ffrozp);
  post.cmc = L.r(L.cmc); post.sneqv = L.r(L.sneqv); post.snowhk = L.r(L.snowhk); post.albedok = L.r(L.albedok);
  post.albbrd = L.r(L.albbrd); post.z0brd = L.r(L.z0brd); post.emissi = L.r(L.emissi); post.z0k = L.r(L.z0k);
  post.t1 = L.r(L.t1); post.sheat = L.r(L.sheat); post.etp = L.r(L.etp); post.eta_kin = L.r(L.eta_kin);
  post.eta = L.r(L.eta); post.ssoil = L.r(L.ssoil); post.sncovr = L.r(L.sncovr); post.snotime1 = L.r(L.snotime1);
  post.q1 = L.r(L.q1); post.flx1 = L.r(L.flx1); post.flx2 = L.r(L.flx2); post.flx3 = L.r(L.flx3);
  post.soilw = L.r(L.soilw); post.soilm = L.r(L.soilm); post.runoff1 = L.r(L.runoff1);
  post.runoff2 = L.r(L.runoff2); post.snomlt = L.r(L.snomlt); post.xlai = L.r(L.xlai);
  post.cqs2 = c2(a.cqs2);
  post.smc = L.r3(L.smc); post.stc = L.r3(L.stc); post.swc = L.r3(L.swc); post.smav = L.r3(L.smav);
  post.lai = w2(a.lai); post.canwat = w2(a.canwat); post.snow = w2(a.snow); post.snowh = w2(a.snowh);
  post.albedo = w2(a.albedo); post.albbck = w2(a.albbck); post.z0 = w2(a.z0); post.emiss = w2(a.emiss);
  post.znt = w2(a.znt); post.tsk = w2(a.tsk); post.hfx = w2(a.hfx); post.potevp = w2(a.potevp);
  post.qfx = w2(a.qfx); post.lh = w2(a.lh); post.grdflx = w2(a.grdflx); post.snowc = w2(a.snowc);
  post.chs2 = w2(a.chs2); post.snotime = w2(a.snotime); post.qsfc = w2(a.qsfc); post.noahres = w2(a.noahres);
  post.smstav = w2(a.smstav); post.smstot = w2(a.smstot); post.sfcrunoff = w2(a.sfcrunoff);
  post.udrunoff = w2(a.udrunoff); post.acsnow = w2(a.acsnow); post.acsnom = w2(a.acsnom);
  post.snopcx = w2(a.snopcx);
  post.smois = v3(a.smois, ni); post.tslb = v3(a.tslb, ni); post.sh2o = v3(a.sh2o, ni);
  post.smcrel = v3(a.smcrel, ni);
  post.dt = a.dt; post.nsoil = nsoil;

  // JLOOP / ILOOP (:744-1593): one thread per point (the points are independent)
  const int ims = a.ims;
  parallel_for(a.its, a.ite, [=] CPPMPAS_LAMBDA(int i) {
    const int ii = i - ims;
    pre(ii, 0);
    nk(ii, 0);   // returns at once unless ICE = 0 on land (active)
    gl(ii, 0);   // returns at once unless ICE = -1
    post(ii, 0); // returns at once on open water and sea ice (CYCLE ILOOP, :1073)
  });
}

// module_sf_noah_seaice_drv.F:18-517 for the points its..ite
void seaice_noah(const SeaiceArgs& a) {
  if (a.seaice_albedo_opt != 0 || a.seaice_thickness_opt != 0 || a.seaice_snowdepth_opt != 0)
    throw Error("module_sf_noah_seaice_drv.F: seaice_albedo_opt, seaice_thickness_opt and seaice_snowdepth_opt "
                "other than 0 are not ported");
  if (a.sf_urban_physics != 0) throw Error("module_sf_noah_seaice_drv.F: sf_urban_physics /= 0 is not ported");
  if (a.kts != a.kms) throw Error("module_sf_noah_seaice_drv.F: kts /= kms is not ported");
  // one layer: SHFLX / SMFLX / the ice and glacial HRT read layer 2 of zsoil, stc, smc (undefined in MPAS)
  if (a.num_soil_layers < 2 || a.num_soil_layers > noah::kMaxSoilLevels)
    throw Error("module_sf_noah_seaice_drv.F: num_soil_layers = " + std::to_string(a.num_soil_layers) +
                " is outside 2.." + std::to_string(noah::kMaxSoilLevels));
  const int ni = a.ime - a.ims + 1;
  auto c2 = [](const real* p) { return v2<const real>(p); };
  auto w2 = [](real* p) { return v2<real>(p); };
  KSeaiceNoah si{};
  si.t3d = v3<const real>(a.t3d, ni); si.qv3d = v3<const real>(a.qv3d, ni); si.p8w3d = v3<const real>(a.p8w3d, ni);
  si.dz8w = v3<const real>(a.dz8w, ni);
  si.glw = c2(a.glw); si.swdown = c2(a.swdown); si.rainbl = c2(a.rainbl); si.snoalb2d = c2(a.snoalb2d);
  si.qgh = c2(a.qgh); si.xice = c2(a.xice); si.rib = c2(a.rib); si.sr = c2(a.sr);
  si.emiss = w2(a.emiss); si.albedo = w2(a.albedo); si.z02d = w2(a.z02d); si.tsk = w2(a.tsk);
  si.snow = w2(a.snow); si.snowc = w2(a.snowc); si.snowh2d = w2(a.snowh2d); si.chs = w2(a.chs);
  si.chs2 = w2(a.chs2); si.cqs2 = w2(a.cqs2); si.znt = w2(a.znt); si.lh = w2(a.lh); si.hfx = w2(a.hfx);
  si.qfx = w2(a.qfx); si.potevp = w2(a.potevp); si.grdflx = w2(a.grdflx); si.qsfc = w2(a.qsfc);
  si.acsnow = w2(a.acsnow); si.acsnom = w2(a.acsnom); si.snopcx = w2(a.snopcx);
  si.sfcrunoff = w2(a.sfcrunoff); si.noahres = w2(a.noahres); si.tslb = v3(a.tslb, ni);
  si.op.albedo_default = a.seaice_albedo_default;
  si.op.snowdepth_max = a.seaice_snowdepth_max;
  si.op.snowdepth_min = a.seaice_snowdepth_min;
  si.thickness_opt = a.seaice_thickness_opt; si.thickness_default = a.seaice_thickness_default;
  si.dt = a.dt; si.xice_threshold = a.xice_threshold; si.nsoil = a.num_soil_layers; si.frpcpn = a.frpcpn;
  // SEAICE_JLOOP / SEAICE_ILOOP (:211-514)
  const int ims = a.ims;
  parallel_for(a.its, a.ite, [=] CPPMPAS_LAMBDA(int i) { si(i - ims, 0); });
}

}  // namespace cppmpas::phys::noahdrv
