# SPDX-License-Identifier: GPL-3.0-or-later
"""Synthetic columns for the surface layer / YSU PBL column check, and the comparison of the two programs' outputs.

  columns.py make <in.bin> <nCellsSolve> <nVertLevels> <nsteps> <seed> <frac_seaice> <ysu_pblmix> <cu> <radt>
  columns.py compare <ref out.bin> <cpp out.bin>

cu: 0 'off', 1 'cu_ntiedtke', 2 'cu_kain_fritsch', 3 'cu_grell_freitas' (the convection branches of
physics_get_tend_work); radt: 0 'off', 1 'rrtmg_lw' / 'rrtmg_sw'.

Both files are native-endian streams of records: name (32 bytes, blank padded), int32 kind (0 float64, 1 int32),
int32 ndim, int32 dims[4], then the data in Fortran order. in.bin starts with the records "hdr" (int32: nCellsSolve,
nCells, nEdgesSolve, nEdges, nVertLevels, nsteps, frac_seaice, ysu_pblmix, cu, radt) and "dtl" (float64: dt,
config_len_disp), then the pool arrays (their Registry names; the harnesses register each in its pool). The outputs
are, in order: the four sf_sfclayrev tables; then per step the snapshots "A" (after MPAS_to_physics: every pool array
and the mpas_atmphys_vars arrays, names starting with "@"), "B" (after driver_sfclayer), "C" (after driver_pbl), "D"
(after physics_get_tend with rk_step = 1, dynamics_substep = 1) and "E" (rk_step = 2), each a marker record then
every pool array.

nCells = nCellsSolve + 2 (two halo cells: their tendencies are read by tend_toEdges, never computed), nEdges =
3 nCells, nEdgesSolve = nEdges - 3. Cell i (1-based) gets regime (i - 1) % 7:
  0 land, daytime convective (hot ground, light wind)       1 land, stable night (cold ground, calm)
  2 sea, unstable (warm sea, moist)                          3 sea, stable (cold sea, strong wind)
  4 sea ice, xice in [0.3, 1] (frac_seaice blend, xice >= 0.5 keeps the regime of the last step)
  5 land, stable air over a cooler ground with mol < 0 from the previous step (br = min(br, 0) = 0: forced
    convection, regime 3), negative and -0.0 mixing ratios from the transport
  6 sea, stratocumulus: cloud water at the top of a well-mixed layer, radiative cooling there (YSU's top-down
    mixing with config_ysu_pblmix), ust < 0.001 in the first step; xice = 0.01 (below the threshold) or 0.02 (on it)
Temperatures, pressures and densities are hydrostatic and consistent; amplitudes are perturbed per cell.
"""
import sys

import numpy as np

RGAS, RV, G = 287.0, 461.6, 9.80616
CP = 7.0 * RGAS / 2.0
P0 = 1.0e5
NAMELEN = 32


def put(o, name, a, kind=0):
    a = np.asarray(a)
    dims = list(a.shape[::-1]) if a.ndim else [1]  # numpy (C order, last fastest) = Fortran reversed
    hdr = name.ljust(NAMELEN).encode()
    o.write(hdr)
    np.array([kind, len(dims)] + dims + [1] * (4 - len(dims)), dtype="<i4").tofile(o)
    np.ascontiguousarray(a, dtype="<f8" if kind == 0 else "<i4").tofile(o)


def records(path):
    out = []
    with open(path, "rb") as f:
        while True:
            h = f.read(NAMELEN)
            if not h:
                break
            name = h.decode().strip()
            meta = np.frombuffer(f.read(24), dtype="<i4")
            kind, nd, dims = int(meta[0]), int(meta[1]), [int(x) for x in meta[2:6]]
            n = int(np.prod(dims[:nd]))
            data = np.frombuffer(f.read(n * (8 if kind == 0 else 4)), dtype="<u8" if kind == 0 else "<i4")
            out.append((name, kind, dims[:nd], data))
    return out


def make(path, nC, nV, nsteps, seed, seaice, pblmix, cu, radt):
    rng = np.random.default_rng(seed)
    nCt = nC + 2
    nE = 3 * nCt
    nEs = nE - 3
    ztop = 25000.0
    eta = np.linspace(0.0, 1.0, nV + 1) ** 1.6  # finer near the ground, ~30 m first layer
    shape_c = (nCt + 1,)
    zg = np.empty((nCt + 1, nV + 1))
    zz = np.empty((nCt + 1, nV))
    names3 = ["rho_zz", "theta_m", "exner", "pressure_base", "pressure_p", "rtheta_base", "rtheta_p",
              "uReconstructZonal", "uReconstructMeridional"]
    f = {n: np.empty((nCt + 1, nV)) for n in names3}
    w = np.empty((nCt + 1, nV + 1))
    sc = np.zeros((nCt + 1, nV, 6))
    sfc = {n: np.empty(shape_c) for n in ["skintemp", "xland", "sst", "xice"]}
    d1 = {n: np.empty(shape_c) for n in
          ["plrad", "hpbl", "mavail", "br", "cpm", "chs", "chs2", "cqs2", "flhc", "flqc", "gz1oz0", "hfx", "qfx",
           "qgh", "qsfc", "lh", "mol", "psih", "psim", "regime", "rmol", "ust", "ustm", "wspd", "znt", "zol", "q2",
           "t2m", "th2m", "u10", "v10", "cd", "cda", "ck", "cka", "fh", "fm", "delta", "wstar"]}
    surface_pressure = np.empty(shape_c)
    rthlw = np.zeros((nCt + 1, nV))
    rthsw = np.zeros((nCt + 1, nV))
    for i in range(nCt + 1):
        reg = i % 7  # array index i = cell i + 1
        land = reg in (0, 1, 5)
        hs = rng.uniform(0.0, 1200.0) if land else 0.0
        zg[i] = hs + eta * (ztop - hs) + 5.0 * np.sin(np.pi * eta) * rng.uniform(-1, 1)
        zm = 0.5 * (zg[i, 1:] + zg[i, :-1])
        zz[i] = (ztop / (ztop - hs)) * (1.0 + 0.003 * rng.uniform(-1, 1, nV))
        z = zm - hs
        tsfc = {0: 300.0, 1: 282.0, 2: 299.0, 3: 284.0, 4: 262.0, 5: 290.0, 6: 289.0}[reg] + rng.uniform(-2, 2)
        lapse = np.full(nV, 0.0065)
        if reg in (0, 2):  # well mixed (dry adiabatic) below the PBL top
            lapse[z < rng.uniform(800, 1800)] = 0.0098
        if reg in (1, 3, 4, 5):  # surface inversion
            lapse[z < rng.uniform(150, 400)] = -0.012 * rng.uniform(0.5, 1.5)
        if reg == 6:  # mixed layer capped by an inversion at ~1 km
            zi = rng.uniform(700, 1100)
            lapse[z < zi] = 0.0098
            lapse[(z >= zi) & (z < zi + 150)] = -0.03
        T = np.empty(nV)
        T[0] = tsfc - lapse[0] * z[0]
        for k in range(1, nV):
            T[k] = T[k - 1] - lapse[k] * (z[k] - z[k - 1])
        T = np.maximum(T, 205.0 + rng.uniform(0, 5))
        ps = rng.uniform(98000.0, 103000.0) * np.exp(-hs / 8000.0)
        p = np.empty(nV)
        pk = ps
        for k in range(nV):
            dz = zm[k] - (hs if k == 0 else zm[k - 1])
            Tk = T[k] if k == 0 else 0.5 * (T[k] + T[k - 1])
            pk = pk * np.exp(-G * dz / (RGAS * Tk))
            p[k] = pk
        es = 610.78 * np.exp(17.27 * (T - 273.16) / (T - 35.86))
        qs = 0.622 * es / (p - es)
        rh = {0: 0.45, 1: 0.8, 2: 0.85, 3: 0.9, 4: 0.85, 5: 0.6, 6: 0.97}[reg] * rng.uniform(0.9, 1.05, nV)
        qv = np.clip(rh * qs, 1e-8, 0.03)
        qv[zm > 12000.0] = 3e-6
        qc = np.zeros(nV)
        qi = np.zeros(nV)
        if reg == 6:  # stratocumulus deck below the inversion
            qc[(z > zi - 300) & (z < zi)] = rng.uniform(2e-4, 6e-4)
        if reg in (2, 3):
            qc[(z > 2000) & (z < 3000)] = 1e-4 * rng.uniform(0, 1)
            qi[(z > 6000) & (z < 8000)] = 2e-5 * rng.uniform(0, 1)
        qr = 1e-5 * rng.uniform(0, 1, nV) * (z < 3000)
        qsn = 1e-5 * rng.uniform(0, 1, nV) * (T < 270)
        qg = np.zeros(nV)
        mix = [qv, qc, qr, qi, qsn, qg]
        for n in range(6):
            sc[i, :, n] = mix[n]
        if reg == 5:  # negative values and -0.0 from the transport
            sc[i, ::5, 1:4] = -1e-9 * rng.uniform(0, 1, (len(range(0, nV, 5)), 3))
            sc[i, 2, 0] = -1e-7
            sc[i, 4, :] = -0.0
        th = T * (P0 / p) ** (RGAS / CP)
        qvp = np.maximum(qv, 0.0)
        thm = th * (1.0 + RV / RGAS * qvp)
        rho = p / (RGAS * T * (1.0 + RV / RGAS * qvp))  # dry density
        pb = p * (1.0 + rng.uniform(-2e-3, 2e-3, nV))
        f["rho_zz"][i] = rho / zz[i]
        f["theta_m"][i] = thm
        f["exner"][i] = (p / P0) ** (RGAS / CP)
        f["pressure_base"][i] = pb
        f["pressure_p"][i] = p - pb
        f["rtheta_base"][i] = f["rho_zz"][i] * th
        f["rtheta_p"][i] = f["rho_zz"][i] * thm - f["rtheta_base"][i]
        spd = {0: 3.0, 1: 1.0, 2: 6.0, 3: 14.0, 4: 8.0, 5: 4.0, 6: 7.0}[reg] * rng.uniform(0.7, 1.3)
        prof = np.log(1.0 + z / 0.1) / np.log(1.0 + 1000.0 / 0.1)
        dirn = rng.uniform(0, 2 * np.pi) + 0.3 * z / 5000.0
        f["uReconstructZonal"][i] = spd * np.minimum(prof, 1.5 + z / 10000.0) * np.cos(dirn)
        f["uReconstructMeridional"][i] = spd * np.minimum(prof, 1.5 + z / 10000.0) * np.sin(dirn)
        w[i] = 0.05 * rng.uniform(-1, 1, nV + 1)
        w[i, 0] = 0.0
        surface_pressure[i] = ps * rng.uniform(0.999, 1.001)
        # surface: the ground or sea temperature relative to the air
        dT = {0: 6.0, 1: -5.0, 2: 1.5, 3: -2.5, 4: -4.0, 5: -3.0, 6: 0.8}[reg] * rng.uniform(0.6, 1.4)
        sfc["skintemp"][i] = T[0] + 0.0098 * z[0] + dT
        sfc["xland"][i] = 1.0 if land else 2.0
        sfc["sst"][i] = sfc["skintemp"][i] if not land else 285.0
        if reg == 4:
            sfc["sst"][i] = 271.0 + rng.uniform(-1, 1)
        sfc["xice"][i] = {4: rng.uniform(0.3, 1.0), 6: [0.01, 0.02][i % 2]}.get(reg, 0.0)
        # the diag_physics state the previous physics step left
        uns = reg in (0, 2, 6)
        d1["hpbl"][i] = rng.uniform(800, 1800) if uns else rng.uniform(50, 300)
        d1["mavail"][i] = 1.0 if not land else rng.uniform(0.1, 0.5)
        d1["br"][i] = -0.5 * rng.uniform(0, 1) if uns else 0.2 * rng.uniform(0, 1)
        d1["cpm"][i] = CP * (1 + 0.8 * qv[0])
        d1["znt"][i] = rng.uniform(0.05, 0.8) if land else rng.uniform(1e-4, 5e-4)
        d1["ust"][i] = rng.uniform(0.1, 0.6) if not (reg == 6) else 0.0005
        d1["ustm"][i] = d1["ust"][i] * rng.uniform(0.9, 1.1)
        d1["mol"][i] = -rng.uniform(0.05, 0.5) if (uns or reg == 5) else rng.uniform(0.01, 0.2)
        d1["hfx"][i] = rng.uniform(50, 350) if uns else -rng.uniform(5, 60)
        d1["qfx"][i] = rng.uniform(2e-5, 2e-4) if uns else rng.uniform(-1e-6, 2e-5)
        d1["qsfc"][i] = qv[0] * rng.uniform(1.0, 1.3) if (land and i % 3) else 0.0
        d1["regime"][i] = 4.0 if uns else 1.0
        d1["wspd"][i] = spd * 0.4
        for n in ["chs", "chs2", "cqs2", "flhc", "flqc", "gz1oz0", "qgh", "lh", "psih", "psim", "rmol", "zol", "q2",
                  "t2m", "th2m", "cd", "cda", "ck", "cka", "delta", "wstar", "plrad"]:
            d1[n][i] = rng.uniform(-1, 1)
        d1["u10"][i] = f["uReconstructZonal"][i, 0] * 0.8
        d1["v10"][i] = f["uReconstructMeridional"][i, 0] * 0.8
        d1["fm"][i] = rng.uniform(4.0, 12.0)
        d1["fh"][i] = d1["fm"][i] * rng.uniform(0.8, 1.2)
        if radt:
            rthlw[i] = -2e-5 * rng.uniform(0.5, 1.5, nV)
            rthsw[i] = 1e-5 * rng.uniform(0, 1, nV)
            if reg == 6:
                rthlw[i][(z > zi - 120) & (z < zi)] = -5e-4 * rng.uniform(0.7, 1.3)
        else:
            rthlw[i] = 1e-5 * rng.uniform(-1, 1, nV)  # read by YSU (rthraten) whatever the radiation option
            rthsw[i] = 1e-5 * rng.uniform(-1, 1, nV)
    # edges
    coe = rng.integers(1, nCt + 1, size=(nE + 1, 2)).astype(np.int32)
    env = rng.normal(size=(nE + 1, 3))
    env /= np.linalg.norm(env, axis=1, keepdims=True)
    east = rng.normal(size=(nCt + 1, 3))
    east /= np.linalg.norm(east, axis=1, keepdims=True)
    north = rng.normal(size=(nCt + 1, 3))
    north -= (north * east).sum(1, keepdims=True) * east
    north /= np.linalg.norm(north, axis=1, keepdims=True)
    with open(path, "wb") as o:
        put(o, "hdr", np.array([nC, nCt, nEs, nE, nV, nsteps, seaice, pblmix, cu, radt], dtype=np.int32), kind=1)
        put(o, "dtl", np.array([600.0, 240000.0]))
        # mesh
        put(o, "latCell", rng.uniform(-1.5, 1.5, nCt + 1))
        put(o, "lonCell", rng.uniform(0, 6.28, nCt + 1))
        put(o, "fzm", rng.uniform(0.4, 0.6, nV))
        put(o, "fzp", rng.uniform(0.4, 0.6, nV))
        put(o, "rdzw", rng.uniform(1e-3, 1e-2, nV))
        put(o, "zgrid", zg)
        put(o, "zz", zz)
        put(o, "meshDensity", rng.uniform(0.2, 1.0, nCt + 1))
        put(o, "areaCell", rng.uniform(4e9, 6e9, nCt + 1))
        put(o, "cellsOnEdge", coe, kind=1)
        put(o, "edgeNormalVectors", env)
        put(o, "east", east)
        put(o, "north", north)
        # state
        put(o, "rho_zz", f["rho_zz"])
        put(o, "theta_m", f["theta_m"])
        put(o, "w", w)
        put(o, "scalars", sc)
        # diag
        put(o, "surface_pressure", surface_pressure)
        for n in ["exner", "pressure_base", "pressure_p", "rtheta_base", "rtheta_p", "uReconstructZonal",
                  "uReconstructMeridional"]:
            put(o, n, f[n])
        put(o, "rho_edge", rng.uniform(0.1, 1.3, (nE + 1, nV)))
        put(o, "tend_u_phys", rng.uniform(-1, 1, (nE + 1, nV)))
        # sfc_input
        for n in ["skintemp", "xland", "sst", "xice"]:
            put(o, n, sfc[n])
        # diag_physics
        for n, a in d1.items():
            put(o, n, a)
        put(o, "kpbl", rng.integers(-5, 5, nCt + 1).astype(np.int32), kind=1)
        for n in ["kzh", "kzm", "kzq", "exch_h"]:
            put(o, n, rng.uniform(-1, 1, (nCt + 1, nV)))
        # tend_physics
        put(o, "rthratenlw", rthlw)
        put(o, "rthratensw", rthsw)
        for n in ["rublten", "rvblten", "rthblten", "rqvblten", "rqcblten", "rqiblten", "rqsblten", "rncblten",
                  "rniblten", "rnifablten", "rnwfablten"]:
            put(o, n, 1e-3 * rng.uniform(-1, 1, (nCt + 1, nV)))
        for n in ["rucuten", "rvcuten", "rthcuten", "rqvcuten", "rqccuten", "rqrcuten", "rqicuten", "rqscuten"]:
            put(o, n, 1e-4 * rng.uniform(-1, 1, (nCt + 1, nV)))
        put(o, "rublten_Edge", rng.uniform(-1, 1, (nE + 1, nV)))
        put(o, "rucuten_Edge", rng.uniform(-1, 1, (nE + 1, nV)))
        # tend
        put(o, "scalars_tend", rng.uniform(-1, 1, (nCt + 1, nV, 6)))
        # the dynamics' arrays of physics_get_tend
        put(o, "tend_ru_physics", rng.uniform(-1, 1, (nE + 1, nV)))
        put(o, "tend_rtheta_physics", rng.uniform(-1, 1, (nCt + 1, nV)))
        put(o, "tend_rho_physics", rng.uniform(-1, 1, (nCt + 1, nV)))


def compare(ref, cpp):
    a, b = records(ref), records(cpp)
    if [(n, k, d) for n, k, d, _ in a] != [(n, k, d) for n, k, d, _ in b]:
        for x, y in zip(a, b):
            if x[:3] != y[:3]:
                print(f"FAIL: record layout differs: ref {x[:3]} cpp {y[:3]}")
                break
        else:
            print(f"FAIL: record counts differ: ref {len(a)} cpp {len(b)}")
        return 1
    bad, total, phase, step = 0, 0, "", 0
    for (name, kind, dims, x), (_, _, _, y) in zip(a, b):
        if name.startswith("#"):
            phase = name[1:]
            continue
        total += x.size
        d = np.nonzero(x != y)[0]
        if d.size:
            bad += 1
            j = int(d[0])
            if kind == 0:
                vx, vy = x[j:j + 1].view("<f8")[0], y[j:j + 1].view("<f8")[0]
            else:
                vx, vy = x[j], y[j]
            idx = np.unravel_index(j, dims[::-1])[::-1]
            if bad <= 25:
                print(f"{phase} {name}: {d.size} of {x.size} differ, first at {tuple(int(t) + 1 for t in idx)}: "
                      f"ref {vx!r} cpp {vy!r}")
    print(f"{total} values in {len(a)} records, {'IDENTICAL' if not bad else str(bad) + ' records differ'}")
    return 1 if bad else 0


if __name__ == "__main__":
    if sys.argv[1] == "make":
        make(sys.argv[2], *map(int, sys.argv[3:11]))
    else:
        sys.exit(compare(sys.argv[2], sys.argv[3]))
