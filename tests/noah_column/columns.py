# SPDX-License-Identifier: GPL-3.0-or-later
"""Synthetic land points for the Noah column check, and the comparison of the two programs' outputs.

  columns.py make <in.bin> <nCells> <nsteps> <dt> <seed> <landuse> <flags>
  columns.py compare <in.bin> <ref out.bin> <cpp out.bin>

landuse: MODIS (MODIFIED_IGBP_MODIS_NOAH: isurban 13, isice 15, water 17) or USGS (isurban 1, isice 24, water 16).
flags: a string of letters: 'r' config_do_restart, 'a' config_sfc_albedo = false, 's' config_sfc_snowalbedo =
false, 'f' config_frac_seaice = false (xice_threshold 0.5).

in.bin (native-endian stream, read by harness.F90 and test_noah_column.cpp): int32 hdr[10] = nCells, nSoilLevels (4),
nsteps, nV (2), isurban, isice, config_sfc_albedo, config_sfc_snowalbedo, config_do_restart, config_frac_seaice;
float64 dt; char[64] mminlu (blank padded); int32 ivgtyp(nC), isltyp(nC); float64 R1(nC, 53), S2(4, nC, 5),
A3(nC, nV+1, 4), SEA(nC, 9) in Fortran order (the names below).
out.bin: the tables (int32 x 66, then float64 x 1146), then R1 and S2 after init_lsm and after every step.

Point regime i % 12 (land use, soil, snow, frozen soil, sea ice, water):
 0 crop / loam, warm, wet, rain, sunny        6 urban (soil 14 on land -> 7 in lsm)
 1 forest / sandy loam, melting snow pack     7 grassland / sand, dry, night, dew (qv > qsat)
 2 tundra / silt loam, frozen soil, cold snow 8 savanna, snow on warm ground, snowfall (t_p < 273.15)
 3 land ice (isice): glacial                   9 bare / clay (BARE category: SHDFAC = 0), drizzle
 4 sea ice, xice in 0.04..1 (both branches     10 wetland / organic, saturated, heavy rain
   of correct_tsk_over_seaice)                 11 open water with sea ice under the threshold (0.01)
 5 open water (xland 2)
The USGS case maps them to USGS categories, with 25 / 26 / 27 (ISURBAN = 1: SHDFAC = 0) on regimes 9 and 6.
Variant (i // 12) % 4 adds: other soil textures (regimes 0 1 2 7 8; every STAS class but 14 and 16 appears in 145
points); SHDMAX = SHDMIN (var 1); on land ice and sea ice rain on a melting surface (1), dew / deposition under a
stable surface layer with a thin pack and no precipitation (2), freezing rain (3); a top layer below SMCDRY
(regime 7, var 1 and 3).
"""
import sys

import numpy as np

R1 = ['acsnom', 'acsnow', 'canwat', 'chs', 'chs2', 'chklowq', 'cpm', 'cqs2', 'glw', 'grdflx', 'gsw', 'hfx', 'lai',
      'lh', 'noahres', 'potevp', 'qfx', 'qgh', 'qsfc', 'br', 'sfc_albedo', 'sfc_emibck', 'sfc_emiss', 'sfcrunoff',
      'smstav', 'smstot', 'snotime', 'snopcx', 'swddif', 'swddir', 'udrunoff', 'z0', 'znt', 'sr', 'rainncv', 'raincv',
      'q2', 't2m', 'th2m', 'shdmin', 'shdmax', 'snoalb', 'sfc_albbck', 'snow', 'snowc', 'snowh', 'tmn', 'skintemp',
      'vegfra', 'xice', 'xland', 'surface_pressure', 'spare']
S2 = ['dzs', 'sh2o', 'smcrel', 'smois', 'tslb']
A3 = ['dz_p', 't_p', 'qv_p', 'pres2_hyd_p']
SEA = ['cpm_sea', 'chs_sea', 'chs2_sea', 'cqs2_sea', 'hfx_sea', 'qfx_sea', 'qgh_sea', 'qsfc_sea', 'lh_sea']
NSOIL, NV = 4, 2
TI, TR = 66, 16 * 50 + 4 + 10 * 30 + 30 + 12  # table dump: ints, reals
# (veg, soil) per regime
MODIS = {0: (12, 6), 1: (1, 3), 2: (18, 4), 3: (15, 16), 4: (15, 16), 5: (17, 14), 6: (13, 14), 7: (10, 1),
         8: (9, 5), 9: (16, 12), 10: (11, 13), 11: (5, 8)}
USGS = {0: (2, 6), 1: (14, 3), 2: (20, 4), 3: (24, 16), 4: (24, 16), 5: (16, 14), 6: (26, 14), 7: (7, 1),
        8: (10, 5), 9: (25, 12), 10: (17, 13), 11: (15, 8)}
# SOILPARM.TBL STAS: (wilting, max) for the soil moisture range
SMC = {1: (.010, .339), 2: (.028, .421), 3: (.047, .434), 4: (.084, .476), 5: (.061, .484), 6: (.066, .439),
       7: (.069, .404), 8: (.120, .464), 9: (.103, .465), 10: (.100, .406), 11: (.126, .468), 12: (.138, .468),
       13: (.066, .439), 14: (0., 1.), 15: (.006, .20), 16: (.028, .421), 17: (.030, .468), 18: (.006, .200),
       19: (.01, .339)}
SOILS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 17, 18, 19]


def qsat(t, p):
    es = 611.2 * np.exp(17.67 * (t - 273.15) / (t - 29.65))
    return 0.622 * es / (p - es)


def make(path, nC, nsteps, dt, seed, landuse, flags):
    rng = np.random.default_rng(seed)
    cats = MODIS if landuse == 'MODIS' else USGS
    isurban, isice = (13, 15) if landuse == 'MODIS' else (1, 24)
    mminlu = 'MODIFIED_IGBP_MODIS_NOAH' if landuse == 'MODIS' else 'USGS'
    r1 = {n: np.zeros(nC) for n in R1}
    s2 = {n: np.zeros((nC, NSOIL)) for n in S2}
    a3 = {n: np.zeros((nC, NV + 1)) for n in A3}
    sea = {n: np.zeros(nC) for n in SEA}
    ivg = np.zeros(nC, dtype='<i4')
    isl = np.zeros(nC, dtype='<i4')
    u = rng.uniform
    for i in range(nC):
        reg = i % 12
        var = (i // 12) % 4
        veg, soil = cats[reg]
        if reg in (0, 1, 2, 7, 8) and var:  # other soil textures of the same regime
            soil = SOILS[(reg + 7 * (i // 12)) % len(SOILS)]
        ivg[i], isl[i] = veg, soil
        tair = {0: 300., 1: 274., 2: 252., 3: 248., 4: 258., 5: 288., 6: 296., 7: 285., 8: 271.5, 9: 293.,
                10: 297., 11: 266.}[reg] + u(-2, 2)
        tsk = tair + {7: -4., 0: 6., 9: 8.}.get(reg, u(-2, 2))
        if reg in (3, 4) and var == 1:  # rain on melting ice / snow: warm air, skin at 0 C
            tair, tsk = 276. + u(0., 2.), 273.15 + u(0., .5)
        if reg in (3, 4) and var == 3:  # freezing rain: air above 0 C, skin below
            tair, tsk = 274. + u(0., 1.), 271. + u(0., 1.)
        ps = u(85000., 102500.)
        r1['surface_pressure'][i] = ps * (1. + u(-1e-4, 1e-4))
        a3['pres2_hyd_p'][i] = [ps, ps - u(400., 800.), ps - u(1200., 1800.)]
        a3['dz_p'][i] = [u(30., 70.), u(80., 120.), 0.]
        a3['t_p'][i] = [tair, tair - u(0.3, 0.8), 0.]
        rh = {7: 1.05, 2: 0.9, 3: 0.85, 10: 0.95}.get(reg, u(0.4, 0.9))
        if reg in (3, 4) and var == 2:  # deposition (dew) on the ice: moist air, cold skin
            rh, tsk = 1.1, tair - 3.
        qv = rh * qsat(tair, ps)
        a3['qv_p'][i] = [qv, 0.95 * qv, 0.]
        r1['qgh'][i] = qsat(tsk, ps)
        r1['glw'][i] = u(180., 420.) if var != 2 else u(150., 250.)
        alb = {1: .45, 2: .6, 3: .75, 4: .55, 8: .5}.get(reg, u(.12, .3))
        r1['sfc_albedo'][i] = alb
        night = reg in (7,)
        r1['gsw'][i] = 0. if night else u(50., 650.) * (1. - alb)
        r1['sfc_albbck'][i] = {3: .7, 4: .65}.get(reg, u(.1, .25))
        r1['sfc_emibck'][i] = u(.9, .97)
        r1['sfc_emiss'][i] = u(.9, .99)
        r1['snoalb'][i] = u(.5, .8)
        r1['chs'][i] = u(.002, .04)
        r1['chs2'][i] = u(.002, .03)
        r1['cqs2'][i] = u(.002, .03) if reg != 9 else 5e-6  # < 1e-5: the other sfc-diagnostics branch
        r1['cpm'][i] = 1004.5 * (1. + 0.86 * qv)
        r1['br'][i] = u(-1.5, 1.) if var != 2 else u(.1, .8)
        r1['z0'][i] = u(.001, .5)
        r1['znt'][i] = r1['z0'][i] * u(.8, 1.2)
        r1['sr'][i] = {8: 1., 2: .9, 1: .4}.get(reg, 0.)
        r1['rainncv'][i] = {0: u(.5, 3.), 8: u(.5, 2.), 9: .05, 10: u(5., 12.), 1: u(0., 1.), 4: u(0., 1.),
                            3: u(0., .5)}.get(reg, 0.)
        if reg in (3, 4, 1) and var == 2:
            r1['rainncv'][i] = 0.  # no precipitation: snow ageing (ALCALC), no new snow
        r1['raincv'][i] = {0: u(0., 2.), 10: u(2., 6.)}.get(reg, 0.)
        r1['lai'][i] = u(.5, 5.)
        r1['canwat'][i] = u(0., .4) if reg in (0, 10, 1) else 0.
        r1['shdmin'][i] = u(1., 30.)
        r1['shdmax'][i] = r1['shdmin'][i] + (u(20., 60.) if var != 1 else 0.)  # var 1: SHDMAX = SHDMIN
        r1['vegfra'][i] = {9: 0., 6: 0.}.get(reg, u(r1['shdmin'][i] - 5., r1['shdmax'][i] + 5.))
        r1['tmn'][i] = u(270., 290.) if reg not in (2, 3) else u(255., 265.)
        for n in ('acsnom', 'acsnow', 'sfcrunoff', 'udrunoff', 'potevp', 'snopcx', 'noahres', 'smstav', 'smstot',
                  'snotime', 'hfx', 'lh', 'qfx', 'grdflx', 'qsfc', 'q2', 't2m', 'th2m', 'swddif', 'swddir'):
            r1[n][i] = u(0., 5.)
        r1['snotime'][i] = u(0., 3.e5)
        sweq = {1: u(20., 80.), 2: u(40., 120.), 3: u(200., 500.), 4: u(5., 40.), 8: u(2., 10.), 11: u(.5, 2.)}.get(reg, 0.)
        if reg in (3, 4) and var >= 2:
            sweq = u(.02, 60.)  # thin pack: SFLX_GLACIAL's SNEQV < 0.1 m, SFLX_SEAICE's < 0.01 m
        r1['snow'][i] = sweq
        dens = u(.15, .35)
        r1['snowh'][i] = sweq * 0.001 / dens if sweq > 0 else 0.
        r1['snowc'][i] = min(1., u(.6, 1.)) if sweq > 0 else 0.
        if reg == 8:
            r1['snowc'][i] = u(.3, .7)
            tsk = 274. + u(0., 2.)  # T1 > 273.14 with snow, SWDOWN > 10
        r1['skintemp'][i] = tsk
        r1['xland'][i] = 2. if reg in (5, 11) else 1.
        r1['xice'][i] = {4: [1., .5, .15, .07, .04, .3][(i // 12) % 6], 11: .01}.get(reg, 0.)
        wlt, mx = SMC[soil]
        dz = [.1, .3, .6, 1.]
        s2['dzs'][i] = dz
        for k in range(NSOIL):
            if reg == 10:
                sm = mx * u(.97, 1.)
            elif reg == 7:
                sm = wlt + (mx - wlt) * u(.02, .15)
                if k == 0 and var % 2:
                    sm = wlt * u(.6, .95)  # below SMCDRY = WLTSMC: DEVAP's FX = 0
            else:
                sm = wlt + (mx - wlt) * u(.3, .9)
            s2['smois'][i, k] = sm
            s2['smcrel'][i, k] = u(0., 1.)
            if reg in (2, 3, 4, 11):
                tl = 262. + 2.5 * k + u(-2, 2)
            elif reg == 1:
                tl = [272.8, 273.0, 273.6, 275.][k] + u(-.3, .3)
            else:
                tl = tair + (r1['tmn'][i] - tair) * (k + .5) / NSOIL + u(-1, 1)
            s2['tslb'][i, k] = tl
            s2['sh2o'][i, k] = sm * u(.5, 1.)
        sea['cpm_sea'][i] = 1004.5 * (1. + 0.86 * qv) * u(.999, 1.001)
        for n in ('chs_sea', 'chs2_sea', 'cqs2_sea'):
            sea[n][i] = u(.002, .03)
        sea['hfx_sea'][i] = u(-20., 200.)
        sea['qfx_sea'][i] = u(0., 1e-4)
        sea['qgh_sea'][i] = qsat(271.4, ps)
        sea['qsfc_sea'][i] = sea['qgh_sea'][i] * u(.95, 1.)
        sea['lh_sea'][i] = u(0., 250.)
    hdr = [nC, NSOIL, nsteps, NV, isurban, isice, 0 if 'a' in flags else 1, 0 if 's' in flags else 1,
           1 if 'r' in flags else 0, 0 if 'f' in flags else 1]
    with open(path, 'wb') as o:
        np.array(hdr, dtype='<i4').tofile(o)
        np.array([dt], dtype='<f8').tofile(o)
        o.write(mminlu.ljust(64).encode())
        ivg.tofile(o)
        isl.tofile(o)
        np.ascontiguousarray(np.stack([r1[n] for n in R1], axis=0), dtype='<f8').tofile(o)  # (nC, 53) Fortran
        np.ascontiguousarray(np.stack([s2[n] for n in S2], axis=0), dtype='<f8').tofile(o)  # (4, nC, 5)
        np.ascontiguousarray(np.stack([a3[n].T for n in A3], axis=0), dtype='<f8').tofile(o)  # (nC, nV+1, 4)
        np.ascontiguousarray(np.stack([sea[n] for n in SEA], axis=0), dtype='<f8').tofile(o)  # (nC, 9)


def compare(inp, ref, cpp):
    hdr = np.fromfile(inp, dtype='<i4', count=10)
    nC, nsoil, nsteps = int(hdr[0]), int(hdr[1]), int(hdr[2])
    a = np.fromfile(ref, dtype='u1')
    b = np.fromfile(cpp, dtype='u1')
    nst = nC * len(R1) + nsoil * nC * len(S2)
    want = 4 * TI + 8 * TR + 8 * nst * (nsteps + 1)
    if a.size != want or b.size != want:
        print(f"FAIL: sizes ref {a.size} cpp {b.size}, expected {want}")
        return 1
    bad = 0
    ti_a, ti_b = a[:4 * TI].view('<i4'), b[:4 * TI].view('<i4')
    tr_a, tr_b = a[4 * TI:4 * TI + 8 * TR].view('<u8'), b[4 * TI:4 * TI + 8 * TR].view('<u8')
    if (ti_a != ti_b).any() or (tr_a != tr_b).any():
        bad += 1
        print(f"tables differ: {np.count_nonzero(ti_a != ti_b)} integers, {np.count_nonzero(tr_a != tr_b)} reals"
              f" (first real at {np.nonzero(tr_a != tr_b)[0][:1]})")
    ra = a[4 * TI + 8 * TR:].view('<u8').reshape(nsteps + 1, nst)
    rb = b[4 * TI + 8 * TR:].view('<u8').reshape(nsteps + 1, nst)
    names = [(n, nC) for n in R1] + [(n, nsoil * nC) for n in S2]
    for s in range(nsteps + 1):
        off = 0
        for n, sz in names:
            d = np.nonzero(ra[s, off:off + sz] != rb[s, off:off + sz])[0]
            if d.size:
                bad += 1
                j = d[0]
                x, y = ra[s, off + j:off + j + 1].view('<f8')[0], rb[s, off + j:off + j + 1].view('<f8')[0]
                print(f"{'init' if s == 0 else 'step ' + str(s)} {n}: {d.size} of {sz} differ, first at {j}: "
                      f"ref {x!r} cpp {y!r}")
            off += sz
    f0 = ra[0].view('<f8')
    fn = ra[-1].view('<f8')
    changed = np.count_nonzero(ra[-1] != ra[0])
    print(f"{nC} points, {nsteps} steps: {ra.size + TI + TR} values, "
          f"{'IDENTICAL' if not bad else str(bad) + ' differences'}; values changed by the steps {changed}"
          f"{'' if np.isfinite(fn).all() and np.isfinite(f0).all() else '; NON-FINITE values'}")
    return 1 if bad else 0


if __name__ == '__main__':
    if sys.argv[1] == 'make':
        make(sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5]), int(sys.argv[6]), sys.argv[7],
             sys.argv[8] if len(sys.argv) > 8 else '-')
    else:
        sys.exit(compare(*sys.argv[2:5]))
