// SPDX-License-Identifier: GPL-3.0-or-later
// See mpas_atmphys_sfclayer.hpp. SL = src/core_atmosphere/physics/mpas_atmphys_driver_sfclayer.F,
// SR = physics/physics_wrf/module_sf_sfclayrev.F, SP = physics/physics_wrf/sf_sfclayrev_pre.F,
// MS = physics/physics_wrf/module_sf_mynn.F, MP = physics/physics_wrf/sf_mynn_pre.F (MPAS v8.4.2).
#include "physics/mpas_atmphys_sfclayer.hpp"

#include <vector>

#include "core/error.hpp"
#include "core/libm.hpp"
#include "physics/mpas_atmphys_constants.hpp"
#include "physics/sf_mynn.hpp"

namespace cppmpas::phys {

namespace {

// The (ims:ime, jms:jme) arrays of one sfclayrev call (SL:953-981 passes the *_p set, SL:983-1012 the *_sea set),
// work arrays a kernel can write (WorkArray), and their kernel view (the same names, as pointers)
#define CPPMPAS_SFCSET(X)                                                                                           \
  X(br) X(cpm) X(chs) X(chs2) X(cqs2) X(flhc) X(flqc) X(gz1oz0) X(hfx) X(qfx) X(qgh) X(qsfc) X(lh) X(mol) X(psim)    \
  X(psih) X(regime) X(rmol) X(ust) X(ustm) X(wspd) X(znt) X(zol) X(q2) X(t2m) X(th2m) X(u10) X(v10) X(cd) X(cda)    \
  X(ck) X(cka) X(fh) X(fm) X(mavail) X(tsk) X(xland)
struct SfcSetView {
#define CPPMPAS_X(a) real* a;
  CPPMPAS_SFCSET(CPPMPAS_X)
#undef CPPMPAS_X
};
struct SfcSet {
#define CPPMPAS_X(a) WorkArray a;
  CPPMPAS_SFCSET(CPPMPAS_X)
#undef CPPMPAS_X
  void allocate(int n) {
#define CPPMPAS_X(a) a.assign(n, real(0));
    CPPMPAS_SFCSET(CPPMPAS_X)
#undef CPPMPAS_X
  }
  SfcSetView view() {
    SfcSetView w{};
#define CPPMPAS_X(a) w.a = a.data();
    CPPMPAS_SFCSET(CPPMPAS_X)
#undef CPPMPAS_X
    return w;
  }
};

template <class T>
sfclayrev::View2<T> v2(WorkArray& a) { return sfclayrev::View2<T>{a.data(), 0}; }
template <class T>
sfclayrev::View2<const T> c2(const WorkArray& a) { return sfclayrev::View2<const T>{a.data(), 0}; }

CPPMPAS_HD inline std::size_t I(int i) { return static_cast<std::size_t>(i - 1); }

// SR:12-279 (sfclayrev) with sf_sfclayrev_pre_run (SP:66-100) and sf_sfclayrev_run (sf_sfclayrev.hpp) for the row
// its..ite: the wrapper copies the arrays to its..ite columns (_hv) and level kts, calls the scheme and copies back;
// isfflx = 1, shalwater_z0 = shalwater_flag = 0, isftcflx = 0, iz0tlnd = 0, scm_force_flux = 0 (SL:30-33,
// mpas_atmphys_vars.F:604); ustm, ck, cka, cd and cda are present.
// tables: the tables where the kernel reads them (KernelCopy)
void sfclayrev_call(const sfclayrev::SFClayRevTables<real>* tables, SfcSet& s, const WorkArray& u1d,
                    const WorkArray& v1d, const WorkArray& t1d, const WorkArray& qv1d, const WorkArray& p1d,
                    const WorkArray& dz1d, const WorkArray& psfc, const WorkArray& dx, const WorkArray& pblh,
                    const WorkArray& lakemask, const WorkArray& water_depth, int its, int ite) {
  using namespace physics_constants;
  sfclayrev::SFClayRevKernel k{};
  k.u1d = c2<real>(u1d);
  k.v1d = c2<real>(v1d);
  k.t1d = c2<real>(t1d);
  k.qv1d = c2<real>(qv1d);
  k.p1d = c2<real>(p1d);
  k.dz8w1d = c2<real>(dz1d);
  k.psfcpa = c2<real>(psfc);
  k.tsk = c2<real>(s.tsk);
  k.xland = c2<real>(s.xland);
  k.dx = c2<real>(dx);
  k.mavail = c2<real>(s.mavail);
  k.pblh = c2<real>(pblh);
  k.lakemask = c2<real>(lakemask);
  k.water_depth = c2<real>(water_depth);
  k.varf = {};
  k.regime = v2<real>(s.regime);
  k.hfx = v2<real>(s.hfx);
  k.qfx = v2<real>(s.qfx);
  k.qsfc = v2<real>(s.qsfc);
  k.mol = v2<real>(s.mol);
  k.rmol = v2<real>(s.rmol);
  k.gz1oz0 = v2<real>(s.gz1oz0);
  k.wspd = v2<real>(s.wspd);
  k.br = v2<real>(s.br);
  k.psim = v2<real>(s.psim);
  k.psih = v2<real>(s.psih);
  k.fm = v2<real>(s.fm);
  k.fh = v2<real>(s.fh);
  k.znt = v2<real>(s.znt);
  k.zol = v2<real>(s.zol);
  k.ust = v2<real>(s.ust);
  k.cpm = v2<real>(s.cpm);
  k.chs2 = v2<real>(s.chs2);
  k.cqs2 = v2<real>(s.cqs2);
  k.chs = v2<real>(s.chs);
  k.flhc = v2<real>(s.flhc);
  k.flqc = v2<real>(s.flqc);
  k.qgh = v2<real>(s.qgh);
  k.ustm = v2<real>(s.ustm);
  k.lh = v2<real>(s.lh);
  k.u10 = v2<real>(s.u10);
  k.v10 = v2<real>(s.v10);
  k.th2 = v2<real>(s.th2m);
  k.t2 = v2<real>(s.t2m);
  k.q2 = v2<real>(s.q2);
  k.ck = v2<real>(s.ck);
  k.cka = v2<real>(s.cka);
  k.cd = v2<real>(s.cd);
  k.cda = v2<real>(s.cda);
  k.tables = tables;
  k.isfflx = true;          // SL:30 isfflx = 1
  k.shalwater_z0 = false;   // shalwater_flag = 0 (mpas_atmphys_vars.F:604)
  k.isftcflx = 0;           // SL:31
  k.iz0tlnd = 0;            // SL:32
  k.scm_force_flux = false; // SL:33
  k.kim_tofd = false;       // not in MMM-physics a4baf7f
  k.tofd_factor = 0;
  k.cp = cp;
  k.g = gravity;
  k.rovcp = rcp;
  k.r = R_d;
  k.xlv = xlv;
  k.svp1 = svp1;
  k.svp2 = svp2;
  k.svp3 = svp3;
  k.svpt0 = svpt0;
  k.ep1 = ep_1;
  k.ep2 = ep_2;
  k.karman = karman;
  k.p1000mb = P0;
  k.nx = ite - its + 1;   // sf_sfclayrev_run's i loops run over its..ite
  k.thcon_lanes = true;   // MPAS -O3: loop 30 vectorised (sf_sfclayrev.hpp)
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) { k(i - its, 0); });  // the views start at column its = 1
}

// MS:19-389 (sfclay_mynn) for the row its..ite, with sf_mynn_pre_run (MP:62-165) and sf_mynn_run (sf_mynn.hpp): the
// wrapper copies the arrays to its..ite columns (_hv) and level kts (and kts + 1 of dz8w, u3d, v3d), calls the scheme
// and copies back; ustm, ck, cka, cd and cda are present (MS:298-310, 373-385), spp_pbl = 0 (f_spp false, MS:231-232).
// qcg (MS:258) is copied but sf_mynn_run never reads it; th3d and pi3d are not read at all. The local wstar and qstar
// (MS:215, 368-369) are dropped. One thread per cell i, which sits at ii = i - its of the n = ite - its + 1 lanes of
// sf_mynn_run's vectorised loop.
void sfclay_mynn_call(const sfmynn::PsiTables* tables, SfcSet& s, WorkArray& ch, const PhysicsVarsView& v,
                      const WorkArray& snowh, const WorkArray& psfc, const WorkArray& dx, const WorkArray& pblh,
                      int itimestep, int its, int ite) {
  using namespace physics_constants;
  const SfcSetView p = s.view();
  real* chp = ch.data();
  const real *snowh_p = snowh.data(), *psfc_p = psfc.data(), *dx_p = dx.data(), *pblh_p = pblh.data();
  const sfmynn::SfMynnConsts c{cp, gravity, rcp, R_d, xlv, svp1, svp2, svp3, svpt0, ep_1, ep_2, karman};  // SL:1016-1044
  const int n = ite - its + 1;
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {
    const std::size_t j = I(i);
    sfmynn::SfMynnColumnIn in;
    in.dz1 = v.dz_p(i, 1);  // MP:127-140 (kts = 1)
    in.u1 = v.u_p(i, 1);
    in.v1 = v.v_p(i, 1);
    in.qv1 = v.qv_p(i, 1);
    in.p1 = v.pres_hyd_p(i, 1);
    in.t1 = v.t_p(i, 1);
    in.rho1 = v.rho_p(i, 1);
    in.dz2 = v.dz_p(i, 2);
    in.u2 = v.u_p(i, 2);
    in.v2 = v.v_p(i, 2);
    in.mavail = p.mavail[j];  // MS:252-261
    in.pblh = pblh_p[j];
    in.xland = p.xland[j];
    in.tsk = p.tsk[j];
    in.psfcpa = psfc_p[j];
    in.snowh = snowh_p[j];
    in.dx = dx_p[j];
    sfmynn::SfMynnColumnInOut o;
    o.ust = p.ust[j];  // MS:240-245
    o.mol = p.mol[j];
    o.qsfc = p.qsfc[j];
    o.regime = p.regime[j];  // MS:264-285
    o.hfx = p.hfx[j];
    o.qfx = p.qfx[j];
    o.lh = p.lh[j];
    o.rmol = p.rmol[j];
    o.qgh = p.qgh[j];
    o.znt = p.znt[j];
    o.zol = p.zol[j];
    o.cpm = p.cpm[j];
    o.chs2 = p.chs2[j];
    o.cqs2 = p.cqs2[j];
    o.chs = p.chs[j];
    o.ch = chp[j];
    o.flhc = p.flhc[j];
    o.flqc = p.flqc[j];
    o.gz1oz0 = p.gz1oz0[j];
    o.wspd = p.wspd[j];
    o.br = p.br[j];
    o.psim = p.psim[j];
    o.psih = p.psih[j];
    o.u10 = o.v10 = o.th2 = o.t2 = o.q2 = 0.;  // MS:288-295
    o.ck = o.cka = o.cd = o.cda = 0.;      // MS:298-305
    o.ustm = p.ustm[j];                     // MS:306-310
    sfmynn::sf_mynn_column(in, o, c, itimestep, i - its, n, *tables);
    p.regime[j] = o.regime;  // MS:335-359
    p.hfx[j] = o.hfx;
    p.qfx[j] = o.qfx;
    p.lh[j] = o.lh;
    p.mol[j] = o.mol;
    p.rmol[j] = o.rmol;
    p.qgh[j] = o.qgh;
    p.qsfc[j] = o.qsfc;
    p.znt[j] = o.znt;
    p.zol[j] = o.zol;
    p.ust[j] = o.ust;
    p.cpm[j] = o.cpm;
    p.chs2[j] = o.chs2;
    p.cqs2[j] = o.cqs2;
    p.chs[j] = o.chs;
    chp[j] = o.ch;
    p.flhc[j] = o.flhc;
    p.flqc[j] = o.flqc;
    p.gz1oz0[j] = o.gz1oz0;
    p.wspd[j] = o.wspd;
    p.br[j] = o.br;
    p.psim[j] = o.psim;
    p.psih[j] = o.psih;
    p.u10[j] = o.u10;  // MS:362-370
    p.v10[j] = o.v10;
    p.th2m[j] = o.th2;
    p.t2m[j] = o.t2;
    p.q2[j] = o.q2;
    p.ck[j] = o.ck;  // MS:373-385
    p.cka[j] = o.cka;
    p.cd[j] = o.cd;
    p.cda[j] = o.cda;
    p.ustm[j] = o.ustm;
  });
}

}  // namespace

void init_sfclayer(const std::string& sfclayer_scheme, sfclayrev::SFClayRevTables<real>& tables) {
  if (sfclayer_scheme == "sf_monin_obukhov_rev") {  // SL:826-827
    sfclayrev::init_sfclayrev_tables_o3(tables);
    return;
  }
  if (sfclayer_scheme == "sf_mynn") {  // SL:829-830: sf_mynn_init -> psi_init(psi_opt = 0)
    sfmynn::init_sfmynn_tables(tables);
    return;
  }
  if (sfclayer_scheme == "off") return;
  throw Error("init_sfclayer: config_sfclayer_scheme = '" + sfclayer_scheme + "' has no verified port");
}

// One kernel per loop of SL (one thread per cell i).
void driver_sfclayer(int itimestep, const sfclayrev::SFClayRevTables<real>& tables, SfclayerPools& mp,
                     const PhysicsVars& vp, int its, int ite, SfclayerSeaArrays* sea_out) {
  using namespace physics_constants;
  const bool rev = mp.sfclayer_scheme == "sf_monin_obukhov_rev", mynn = mp.sfclayer_scheme == "sf_mynn";
  if (!rev && !mynn)
    throw Error("driver_sfclayer: config_sfclayer_scheme = '" + mp.sfclayer_scheme + "' has no verified port");
  if (its != 1) throw Error("driver_sfclayer: its /= 1");
  if (vp.ime < ite) throw Error("driver_sfclayer: MPAS_to_physics has not filled its..ite");
  // sf_mynn_pre reads level kts+1 of dz8w, u3d and v3d (sf_mynn_pre.F:127-140): at least 2 levels
  if (mynn && vp.kme < 3) throw Error("driver_sfclayer: sf_mynn with fewer than 2 levels (it reads level 2)");
  if (mynn && (!mp.qcg.p || !mp.snowh.p || !mp.ch.p))
    throw Error("driver_sfclayer: sf_mynn's qcg, snowh and ch (SL:568-587) are not set");
  if (mp.config_frac_seaice && (!mp.sst.p || !mp.xice.p))
    throw Error("driver_sfclayer: config_frac_seaice needs sst and xice (SL:352-592)");
  const SfclayerArrays m = mp;
  const PhysicsVarsView v = vp.view();
  const bool seaice = m.config_frac_seaice;
  const int n = ite;  // allocate_sfclayer: (ims:ime, jms:jme), ims = 1
  // allocate_sfclayer (SL:97-227)
  SfcSet ps, seas;
  ps.allocate(n);
  if (seaice) seas.allocate(n);
  WorkArray dx_a(n), hpbl_a(n), sst_a, xice_a, regime_hold_a, waterdepth_a, lakedepth_a, lakemask_a;
  WorkArray snowh_a, ch_a, qcg_a, ch_sea_a;  // SL:215-221 (sf_mynn)
  if (rev) {  // SL:206-209
    waterdepth_a.assign(n, 0);
    lakedepth_a.assign(n, 0);
    lakemask_a.assign(n, 0);
  }
  if (mynn) {
    snowh_a.assign(n, 0);
    ch_a.assign(n, 0);
    qcg_a.assign(n, 0);
    if (seaice) ch_sea_a.assign(n, 0);
  }
  if (seaice) {
    sst_a.assign(n, 0);
    xice_a.assign(n, 0);
    regime_hold_a.assign(n, 0);
  }
  const SfcSetView p = ps.view(), sea = seas.view();
  real *dx_p = dx_a.data(), *hpbl_p = hpbl_a.data(), *sst_p = sst_a.data(), *xice_p = xice_a.data(),
       *regime_hold = regime_hold_a.data(), *waterdepth_p = waterdepth_a.data(), *lakedepth_p = lakedepth_a.data(),
       *lakemask_p = lakemask_a.data();
  real *snowh_p = snowh_a.data(), *ch_p = ch_a.data(), *qcg_p = qcg_a.data(), *ch_sea = ch_sea_a.data();

  // sfclayer_from_MPAS (SL:352-592)
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // SL:422-470
    dx_p[I(i)] = m.config_len_disp / libm::pow(m.meshDensity(i), 0.25);
    hpbl_p[I(i)] = m.hpbl(i);
    p.mavail[I(i)] = m.mavail(i);
    p.tsk[I(i)] = m.skintemp(i);
    p.xland[I(i)] = m.xland(i);
    p.br[I(i)] = m.br(i);
    p.cpm[I(i)] = m.cpm(i);
    p.chs[I(i)] = m.chs(i);
    p.chs2[I(i)] = m.chs2(i);
    p.cqs2[I(i)] = m.cqs2(i);
    p.flhc[I(i)] = m.flhc(i);
    p.flqc[I(i)] = m.flqc(i);
    p.gz1oz0[I(i)] = m.gz1oz0(i);
    p.hfx[I(i)] = m.hfx(i);
    p.qfx[I(i)] = m.qfx(i);
    p.qgh[I(i)] = m.qgh(i);
    p.qsfc[I(i)] = m.qsfc(i);
    p.lh[I(i)] = m.lh(i);
    p.mol[I(i)] = m.mol(i);
    p.psim[I(i)] = m.psim(i);
    p.psih[I(i)] = m.psih(i);
    p.regime[I(i)] = m.regime(i);
    p.rmol[I(i)] = m.rmol(i);
    p.ust[I(i)] = m.ust(i);
    p.wspd[I(i)] = m.wspd(i);
    p.znt[I(i)] = m.znt(i);
    p.zol[I(i)] = m.zol(i);
    p.q2[I(i)] = 0.;
    p.t2m[I(i)] = 0.;
    p.th2m[I(i)] = 0.;
    p.u10[I(i)] = 0.;
    p.v10[I(i)] = 0.;
    p.cd[I(i)] = 0.;
    p.cda[I(i)] = 0.;
    p.ck[I(i)] = 0.;
    p.cka[I(i)] = 0.;
    p.ustm[I(i)] = m.ustm(i);
  });
  if (seaice)  // SL:471-534
    parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {
      sst_p[I(i)] = m.sst(i);
      xice_p[I(i)] = m.xice(i);
      sea.mavail[I(i)] = m.mavail(i);
      sea.tsk[I(i)] = m.skintemp(i);
      sea.xland[I(i)] = m.xland(i);
      sea.br[I(i)] = m.br(i);
      sea.cpm[I(i)] = m.cpm(i);
      sea.chs[I(i)] = m.chs(i);
      sea.chs2[I(i)] = m.chs2(i);
      sea.cqs2[I(i)] = m.cqs2(i);
      sea.flhc[I(i)] = m.flhc(i);
      sea.flqc[I(i)] = m.flqc(i);
      sea.gz1oz0[I(i)] = m.gz1oz0(i);
      sea.lh[I(i)] = m.lh(i);
      sea.hfx[I(i)] = m.hfx(i);
      sea.qfx[I(i)] = m.qfx(i);
      sea.mol[I(i)] = m.mol(i);
      sea.psim[I(i)] = m.psim(i);
      sea.psih[I(i)] = m.psih(i);
      sea.qgh[I(i)] = m.qgh(i);
      sea.rmol[I(i)] = m.rmol(i);
      sea.regime[I(i)] = m.regime(i);
      sea.ust[I(i)] = m.ust(i);
      sea.ustm[I(i)] = m.ustm(i);
      sea.wspd[I(i)] = m.wspd(i);
      sea.zol[I(i)] = m.zol(i);
      sea.znt[I(i)] = m.znt(i);
      regime_hold[I(i)] = m.regime(i);
      sea.cd[I(i)] = 0.;
      sea.cda[I(i)] = 0.;
      sea.ck[I(i)] = 0.;
      sea.cka[I(i)] = 0.;
      sea.qsfc[I(i)] = 0.;
      sea.q2[I(i)] = 0.;
      sea.t2m[I(i)] = 0.;
      sea.th2m[I(i)] = 0.;
      sea.u10[I(i)] = 0.;
      sea.v10[I(i)] = 0.;
      if (xice_p[I(i)] >= m.xice_threshold && xice_p[I(i)] <= 1.) {  // SL:520-530
        sea.xland[I(i)] = 2.;
        sea.mavail[I(i)] = 1.;
        sea.znt[I(i)] = 0.0001;
        sea.tsk[I(i)] = sfclayrev::c_max(sst_p[I(i)], real(271.4));
      } else {
        sea.xland[I(i)] = p.xland[I(i)];
        sea.mavail[I(i)] = p.mavail[I(i)];
        sea.znt[I(i)] = p.znt[I(i)];
        sea.tsk[I(i)] = p.tsk[I(i)];
      }
    });
  if (rev) {
    parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // SL:541-550 (sf_monin_obukhov_rev)
      p.fh[I(i)] = m.fh(i);
      p.fm[I(i)] = m.fm(i);
      if (seaice) {
        sea.fh[I(i)] = m.fh(i);
        sea.fm[I(i)] = m.fm(i);
      }
    });
    parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // SL:556-562
      waterdepth_p[I(i)] = 0.;
      lakedepth_p[I(i)] = 0.;
      lakemask_p[I(i)] = 0.;
    });
  }
  if (mynn)
    parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // SL:568-587 (sf_mynn)
      snowh_p[I(i)] = m.snowh(i);
      qcg_p[I(i)] = m.qcg(i);
      ch_p[I(i)] = m.ch(i);
      if (seaice) ch_sea[I(i)] = m.ch(i);
    });

  // SL:883: dx = sqrt(maxval(areaCell)) is an argument of neither scheme (both take dx_p)
  int initflag = 1;  // SL:885-886
  if (m.config_do_restart || itimestep > 1) initflag = 0;

  if (mynn) {  // SL:1014-1077: sfclay_mynn on the *_p set, then on the *_sea set (snowh_p, qcg_p, hpbl_p shared)
    KernelCopy<sfmynn::PsiTables> tables_k;
    const sfmynn::PsiTables* tk = kernel_ptr(tables, tables_k);
    sfclay_mynn_call(tk, ps, ch_a, v, snowh_a, vp.psfc_p.a, dx_a, hpbl_a, initflag, its, ite);
    if (seaice) sfclay_mynn_call(tk, seas, ch_sea_a, v, snowh_a, vp.psfc_p.a, dx_a, hpbl_a, initflag, its, ite);
  }

  if (rev) {  // SL:951-1013: sfclayrev; the inputs of level kts (SR:156-165, SP:89-96)
    WorkArray u1d_a(n), v1d_a(n), t1d_a(n), qv1d_a(n), p1d_a(n), dz1d_a(n), psfc_a(n);
    {
      real *u1d = u1d_a.data(), *v1d = v1d_a.data(), *t1d = t1d_a.data(), *qv1d = qv1d_a.data(), *p1d = p1d_a.data(),
           *dz1d = dz1d_a.data(), *psfc = psfc_a.data();
      parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {
        u1d[I(i)] = v.u_p(i, 1);
        v1d[I(i)] = v.v_p(i, 1);
        t1d[I(i)] = v.t_p(i, 1);
        qv1d[I(i)] = v.qv_p(i, 1);
        p1d[I(i)] = v.pres_hyd_p(i, 1);
        dz1d[I(i)] = v.dz_p(i, 1);
        psfc[I(i)] = v.psfc_p(i);
      });
    }
    KernelCopy<sfclayrev::SFClayRevTables<real>> tables_k;
    const sfclayrev::SFClayRevTables<real>* tk = kernel_ptr(tables, tables_k);
    sfclayrev_call(tk, ps, u1d_a, v1d_a, t1d_a, qv1d_a, p1d_a, dz1d_a, psfc_a, dx_a, hpbl_a, lakemask_a, waterdepth_a,
                   its, ite);
    if (seaice)
      sfclayrev_call(tk, seas, u1d_a, v1d_a, t1d_a, qv1d_a, p1d_a, dz1d_a, psfc_a, dx_a, hpbl_a, lakemask_a,
                     waterdepth_a, its, ite);
  }

  // sfclayer_to_MPAS (SL:595-796)
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // SL:675-714
    m.br(i) = p.br[I(i)];
    m.cpm(i) = p.cpm[I(i)];
    m.chs(i) = p.chs[I(i)];
    m.chs2(i) = p.chs2[I(i)];
    m.cqs2(i) = p.cqs2[I(i)];
    m.flhc(i) = p.flhc[I(i)];
    m.flqc(i) = p.flqc[I(i)];
    m.gz1oz0(i) = p.gz1oz0[I(i)];
    m.hfx(i) = p.hfx[I(i)];
    m.lh(i) = p.lh[I(i)];
    m.mol(i) = p.mol[I(i)];
    m.qfx(i) = p.qfx[I(i)];
    m.qgh(i) = p.qgh[I(i)];
    m.qsfc(i) = p.qsfc[I(i)];
    m.psim(i) = p.psim[I(i)];
    m.psih(i) = p.psih[I(i)];
    m.regime(i) = p.regime[I(i)];
    m.rmol(i) = p.rmol[I(i)];
    m.ust(i) = p.ust[I(i)];
    m.wspd(i) = p.wspd[I(i)];
    m.zol(i) = p.zol[I(i)];
    m.znt(i) = p.znt[I(i)];
    m.q2(i) = p.q2[I(i)];
    m.t2m(i) = p.t2m[I(i)];
    m.th2m(i) = p.th2m[I(i)];
    m.u10(i) = p.u10[I(i)];
    m.v10(i) = p.v10[I(i)];
    m.cd(i) = p.cd[I(i)];
    m.cda(i) = p.cda[I(i)];
    m.ck(i) = p.ck[I(i)];
    m.cka(i) = p.cka[I(i)];
    m.ustm(i) = p.ustm[I(i)];
  });
  if (seaice)  // SL:715-748
    parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {
      const real x = m.xice(i);
      if (x >= m.xice_threshold && x <= 1.) {
        auto mix = [&](const real* a, const real* b) { return a[I(i)] * x + (1. - x) * b[I(i)]; };
        m.br(i) = mix(p.br, sea.br);
        m.flhc(i) = mix(p.flhc, sea.flhc);
        m.flqc(i) = mix(p.flqc, sea.flqc);
        m.gz1oz0(i) = mix(p.gz1oz0, sea.gz1oz0);
        m.mol(i) = mix(p.mol, sea.mol);
        m.psih(i) = mix(p.psih, sea.psih);
        m.psim(i) = mix(p.psim, sea.psim);
        m.rmol(i) = mix(p.rmol, sea.rmol);
        m.ust(i) = mix(p.ust, sea.ust);
        m.wspd(i) = mix(p.wspd, sea.wspd);
        m.zol(i) = mix(p.zol, sea.zol);
        if (x >= 0.5) m.regime(i) = regime_hold[I(i)];  // SL:731
        m.q2(i) = mix(p.q2, sea.q2);
        m.t2m(i) = mix(p.t2m, sea.t2m);
        m.th2m(i) = mix(p.th2m, sea.th2m);
        m.u10(i) = mix(p.u10, sea.u10);
        m.v10(i) = mix(p.v10, sea.v10);
        m.cd(i) = mix(p.cd, sea.cd);
        m.cda(i) = mix(p.cda, sea.cda);
        m.ck(i) = mix(p.ck, sea.ck);
        m.cka(i) = mix(p.cka, sea.cka);
        m.ustm(i) = mix(p.ustm, sea.ustm);
      }
    });
  if (rev) {
    parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // SL:755-760
      m.fh(i) = p.fh[I(i)];
      m.fm(i) = p.fm[I(i)];
    });
    if (seaice)  // SL:761-772
      parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {
        const real x = m.xice(i);
        if (x >= m.xice_threshold && x <= 1.) {
          m.fh(i) = p.fh[I(i)] * x + (1. - x) * sea.fh[I(i)];
          m.fm(i) = p.fm[I(i)] * x + (1. - x) * sea.fm[I(i)];
        }
      });
  }
  if (mynn) {
    parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) { m.ch(i) = ch_p[I(i)]; });  // SL:776-780
    if (seaice)  // SL:781-790
      parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {
        const real x = m.xice(i);
        if (x >= m.xice_threshold && x <= 1.) m.ch(i) = ch_p[I(i)] * x + (1. - x) * ch_sea[I(i)];
      });
  }
  if (seaice && sea_out) {  // the module arrays the sea-ice driver reads next
    sea_out->chs = seas.chs;
    sea_out->chs2 = seas.chs2;
    sea_out->cqs2 = seas.cqs2;
    sea_out->cpm = seas.cpm;
    sea_out->hfx = seas.hfx;
    sea_out->qfx = seas.qfx;
    sea_out->qgh = seas.qgh;
    sea_out->qsfc = seas.qsfc;
    sea_out->lh = seas.lh;
  }
}

}  // namespace cppmpas::phys
