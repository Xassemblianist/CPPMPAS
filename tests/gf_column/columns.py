# SPDX-License-Identifier: GPL-3.0-or-later
"""Synthetic columns for the Grell-Freitas column check, and the comparison of the two programs' outputs.

  columns.py make <in.bin> <nCellsSolve> <nVertLevels> <nsteps> <seed> <len_disp> <acrain> <closure_deep>
             <closure_shallow> <dt>
  columns.py compare <ref out.bin> <cpp out.bin>
  columns.py stats <ref out.bin> <nCellsSolve>

len_disp: config_len_disp (dx_p = len_disp / meshDensity**0.25; 240000 is RG's; 3000 makes Grell-Freitas' scale
awareness act: frh = 0.7 and sig = 0.09, so the 1200 code of saturated ascent is reached); acrain: 1 applies the
rainc bucket (l_acrain, config_bucket_rainc = 2 mm); closure_deep / closure_shallow: config_gfconv_closure_deep /
_shallow (Registry defaults 0 and 8); dt: config_dt.

Both files are native-endian streams of records: name (32 bytes, blank padded), int32 kind (0 float64, 1 int32),
int32 ndim, int32 dims[4], then the data in Fortran order. in.bin starts with the records "hdr" (int32: nCellsSolve,
nCells, nEdgesSolve, nEdges, nVertLevels, nsteps, acrain, restart = 0, closure_deep, closure_shallow) and "dtl"
(float64: dt, config_len_disp, config_bucket_rainc), then the pool arrays (their Registry names; the harnesses
register each in its pool). The outputs are, per step, the snapshots "A" (after MPAS_to_physics: every pool array),
"C" (after update_convection_step1, driver_convection and update_convection_step2), "D" (after physics_get_tend with
rk_step = 1, dynamics_substep = 1) and "E" (rk_step = 2), each a marker record then every pool array.

nCells = nCellsSolve + 2 (two halo cells), nEdges = 3 nCells, nEdgesSolve = nEdges - 3. Cell i (1-based) gets regime
(i - 1) % 10:
  0 land, deep convective (hot, moist boundary layer, large CAPE, upward motion, moisture convergence), kpbl 4-12
  1 sea, deep convective (tropical), kpbl 3-9
  2 land, stable and dry (no convection), kpbl 1-3
  3 sea, stable and dry, kpbl 1-4
  4 land, moist mid troposphere with large-scale ascent over a stable boundary layer, cold top (ice), kpbl 2-6
  5 sea, trade-wind cumulus under an inversion at 0.9-1.4 km (shallow convection), kpbl 6-10 (the shallow scheme's
    k22 = kpbl branch); negative and -0.0 mixing ratios from the transport
  6 land, shallow cumulus over a heated dry boundary layer (large surface fluxes, kpbl 2-5), cold column top
  7 as 1, saturated below 3 km with strong ascent (1-4 m/s): with len_disp = 3000 the 1200 code (cup_forcing_ens_3d
    clears it), and large closures (neg_check's limits); negative water vapour at a few levels of 2.5-6 km
  (5 also gets, in half of its cells, a nearly dry level, 1e-9, inside the cloud layer)
  8 as 1, with forcing that stabilises the column (aa1 < aa0: the closures of xff0 < 0 over water)
  9 as 4 over water, with a cold, dry layer at the ground (an elevated cloud; the downdraft-origin search)
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
    o.write(name.ljust(NAMELEN).encode())
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


def make(path, nC, nV, nsteps, seed, len_disp, acrain, cdeep, cshal, dt):
    rng = np.random.default_rng(seed)
    nCt = nC + 2
    nE = 3 * nCt
    nEs = nE - 3
    ztop = 30000.0
    eta = np.linspace(0.0, 1.0, nV + 1) ** 1.5
    sh = (nCt + 1,)
    zg = np.empty((nCt + 1, nV + 1))
    zz = np.empty((nCt + 1, nV))
    names3 = ["rho_zz", "theta_m", "exner", "pressure_base", "pressure_p", "uReconstructZonal",
              "uReconstructMeridional"]
    f = {n: np.empty((nCt + 1, nV)) for n in names3}
    w = np.empty((nCt + 1, nV + 1))
    sc = np.zeros((nCt + 1, nV, 6))
    kpbl = np.ones(sh, dtype=np.int32)
    xland = np.empty(sh)
    hfx = np.empty(sh)
    qfx = np.empty(sh)
    sp = np.empty(sh)
    t3 = {n: np.zeros((nCt + 1, nV)) for n in ["rthdynten", "rqvdynten", "rthblten", "rqvblten", "rublten",
                                                 "rvblten", "rthratenlw", "rthratensw"]}
    for i in range(nCt + 1):
        ext = i % 10  # array index i = cell i + 1
        reg = {7: 1, 8: 1, 9: 4}.get(ext, ext)  # 7, 8, 9: variants of the tropical and the elevated columns
        land = reg in (0, 2, 4, 6) and ext != 9
        hs = rng.uniform(100.0, 1500.0) if land else 0.0
        zg[i] = hs + eta * (ztop - hs) + 3.0 * np.sin(np.pi * eta) * rng.uniform(-1, 1)
        zm = 0.5 * (zg[i, 1:] + zg[i, :-1])
        zz[i] = (ztop / (ztop - hs)) * (1.0 + 0.003 * rng.uniform(-1, 1, nV))
        z = zm - hs
        tsfc = {0: 303.0, 1: 301.5, 2: 278.0, 3: 286.0, 4: 288.0, 5: 297.0, 6: 299.0}[reg] + rng.uniform(-1.5, 1.5)
        lapse = np.full(nV, 0.0065)
        ztp = rng.uniform(14000.0, 16500.0) if reg in (0, 1, 5) else rng.uniform(10000.0, 12000.0)
        if reg in (0, 1):  # deep convective: mixed below ~1 km, moist adiabat-like above
            lapse[z < rng.uniform(600, 1200)] = 0.0085
            lapse[(z > 1200) & (z < ztp)] = 0.0058
        if reg in (2, 3):  # stable: inversion near the ground, stable aloft
            lapse[z < rng.uniform(200, 600)] = -0.008
            lapse[(z > 600) & (z < ztp)] = 0.0045
        if reg == 4:  # cool, stable boundary layer, moist neutral mid troposphere
            lapse[z < 800] = 0.002
            lapse[(z > 2500) & (z < 7000)] = 0.0072
        if reg == 5:  # trade-wind boundary layer: dry adiabatic to cloud base, conditionally unstable to the inversion
            zi = rng.uniform(900, 1400)  # cloud top within 200 hPa of the ground: shallow (cutypen)
            zb = rng.uniform(250, 450)
            lapse[z < zi] = 0.0058
            lapse[z < zb] = 0.0085
            lapse[(z >= zi) & (z < zi + 700)] = -0.02
        if reg == 6:  # heated dry boundary layer, moist layer above it, then stable
            lapse[z < 1500] = 0.0095
            lapse[(z >= 1500) & (z < 3000)] = 0.0062
            lapse[(z >= 3000) & (z < 3600)] = -0.004
        lapse[z >= ztp] = -0.002
        T = np.empty(nV)
        T[0] = tsfc - lapse[0] * z[0]
        for k in range(1, nV):
            T[k] = T[k - 1] - lapse[k] * (z[k] - z[k - 1])
        if ext == 9:
            T[z < 1500] = T[z < 1500] - rng.uniform(3.0, 12.0) * (1.0 - z[z < 1500] / 1500.0)
        T = np.maximum(T, 195.0 + rng.uniform(0, 5))
        ps = rng.uniform(99500.0, 101500.0) * np.exp(-hs / 8000.0)
        p = np.empty(nV)
        pk = ps
        for k in range(nV):
            dz = zm[k] - (hs if k == 0 else zm[k - 1])
            Tk = T[k] if k == 0 else 0.5 * (T[k] + T[k - 1])
            pk = pk * np.exp(-G * dz / (RGAS * Tk))
            p[k] = pk
        es = 610.78 * np.exp(17.27 * (T - 273.16) / (T - 35.86))
        qs = 0.622 * es / (p - es)
        if reg in (0, 1):
            rh = np.where(z < 1500, 0.92, np.where(z < 8000, 0.75, 0.4))
        elif reg in (2, 3):
            rh = np.where(z < 500, 0.6, 0.15)
        elif reg == 4:
            rh = np.where(z < 1500, 0.7, np.where(z < 8000, 0.97, 0.5))
        elif reg == 5:
            rh = np.where(z < zi, 0.92, 0.2)
        else:
            rh = np.where(z < 1200, 0.45, np.where(z < 3000, 0.85, 0.25))
        rh = rh * rng.uniform(0.93, 1.04, nV)
        if ext == 7:  # saturated below 3 km (q_cup / qes_cup >= 0.99 at k22)
            rh = np.where(z < 3000, 1.0, rh)
        if ext == 9:  # elevated: a cold, dry layer at the ground under the moist one
            rh = np.where(z < 1500, 0.3, rh)
        qv = np.clip(rh * qs, 1e-8, 0.03)
        qv[zm > ztp] = 3e-6
        qc = np.zeros(nV)
        qi = np.zeros(nV)
        if reg in (0, 1, 4):
            qc[(z > 1500) & (z < 5000) & (T > 260)] = 1e-4 * rng.uniform(0, 1)
            qi[(T < 250) & (zm < ztp)] = 3e-5 * rng.uniform(0, 1)
        if reg == 5:
            qc[(z > zi - 400) & (z < zi)] = rng.uniform(1e-4, 4e-4)
        qr = 1e-5 * rng.uniform(0, 1, nV) * (z < 3000)
        qsn = 1e-5 * rng.uniform(0, 1, nV) * (T < 270)
        mix = [qv, qc, qr, qi, qsn, np.zeros(nV)]
        for n in range(6):
            sc[i, :, n] = mix[n]
        if reg == 5:  # negative values and -0.0 from the transport, above the trade inversion
            sc[i, 15::5, 1:4] = -1e-9 * rng.uniform(0, 1, (len(range(15, nV, 5)), 3))
            sc[i, 2, 2:4] = -1e-9
            sc[i, 18, 0] = -1e-7
            sc[i, 21, :] = -0.0
            kc = np.nonzero((z > zb) & (z < zi))[0]  # a dry level inside the cloud layer (the moisture limit)
            if kc.size and rng.uniform() < 0.5:
                sc[i, kc[kc.size // 2], 0] = 1e-9
        if ext == 7:  # negative water vapour inside the deep cloud (q2d = 1e-8; neg_check's moisture limit)
            kc = np.nonzero((z > 2500) & (z < 6000))[0][::3]
            sc[i, kc, 0] = -1e-7 * rng.uniform(0, 1, kc.size)
        th = T * (P0 / p) ** (RGAS / CP)
        qvp = np.maximum(sc[i, :, 0], 0.0)
        thm = th * (1.0 + RV / RGAS * qvp)
        rho = p / (RGAS * T * (1.0 + RV / RGAS * qvp))
        pb = p * (1.0 + rng.uniform(-2e-3, 2e-3, nV))
        f["rho_zz"][i] = rho / zz[i]
        f["theta_m"][i] = thm
        f["exner"][i] = (p / P0) ** (RGAS / CP)
        f["pressure_base"][i] = pb
        f["pressure_p"][i] = p - pb
        spd = {0: 4.0, 1: 6.0, 2: 16.0, 3: 9.0, 4: 10.0, 5: 8.0, 6: 5.0}[reg] * rng.uniform(0.7, 1.3)
        prof = np.log(1.0 + z / 0.1) / np.log(1.0 + 1000.0 / 0.1)
        dirn = rng.uniform(0, 2 * np.pi) + 0.4 * z / 6000.0
        f["uReconstructZonal"][i] = spd * np.minimum(prof, 1.2 + z / 8000.0) * np.cos(dirn)
        f["uReconstructMeridional"][i] = spd * np.minimum(prof, 1.2 + z / 8000.0) * np.sin(dirn)
        zw = zg[i] - hs
        if reg in (0, 1):
            w[i] = 0.4 * np.sin(np.pi * np.minimum(zw / ztp, 1.0)) * rng.uniform(0.3, 1.0)
        elif reg == 4:
            w[i] = 0.15 * np.sin(np.pi * np.minimum(zw / 10000.0, 1.0)) * rng.uniform(0.5, 1.0)
        else:
            w[i] = 0.03 * rng.uniform(-1, 1, nV + 1)
        if ext == 7:  # strong ascent (large omega and moisture-convergence closures)
            w[i] = rng.uniform(1.0, 4.0) * np.sin(np.pi * np.minimum(zw / ztp, 1.0))
        w[i, 0] = 0.0
        sp[i] = ps * rng.uniform(0.999, 1.001)
        xland[i] = 1.0 if land else 2.0
        hfx[i] = {0: rng.uniform(150, 400), 1: rng.uniform(5, 30), 2: -rng.uniform(5, 50), 3: rng.uniform(-10, 5),
                  4: rng.uniform(-20, 40), 5: rng.uniform(100, 250), 6: rng.uniform(250, 450)}[reg]
        qfx[i] = {0: rng.uniform(5e-5, 2e-4), 1: rng.uniform(5e-5, 1.5e-4), 2: rng.uniform(0, 1e-5),
                  3: rng.uniform(0, 2e-5), 4: rng.uniform(0, 3e-5), 5: rng.uniform(8e-5, 2e-4),
                  6: rng.uniform(1e-4, 2.5e-4)}[reg]
        kpbl[i] = {0: rng.integers(4, 13), 1: rng.integers(3, 10), 2: rng.integers(1, 4), 3: rng.integers(1, 5),
                   4: rng.integers(2, 7), 5: rng.integers(6, 11), 6: rng.integers(2, 6)}[reg]
        # tendencies from the dynamics, the PBL and the radiation
        conv = reg in (0, 1, 4, 5, 6)
        lowz = np.exp(-z / 1500.0)
        t3["rqvdynten"][i] = (4e-7 * lowz * rng.uniform(0.5, 1.5) if conv else 0.0) + 1e-8 * rng.uniform(-1, 1, nV)
        t3["rthdynten"][i] = -5e-6 * lowz * rng.uniform(0, 1) + 1e-6 * rng.uniform(-1, 1, nV)
        hb = np.exp(-z / 800.0)
        t3["rthblten"][i] = hfx[i] / 300.0 * 1e-3 * hb * rng.uniform(0.5, 1.5)
        t3["rqvblten"][i] = qfx[i] * 1e-3 * hb * rng.uniform(0.5, 1.5) * 10.0
        t3["rublten"][i] = -1e-3 * hb * f["uReconstructZonal"][i] + 1e-5 * rng.uniform(-1, 1, nV)
        t3["rvblten"][i] = -1e-3 * hb * f["uReconstructMeridional"][i] + 1e-5 * rng.uniform(-1, 1, nV)
        if ext == 8:  # forcing that stabilises: drying at low levels, warming aloft (aa1 < aa0)
            t3["rqvdynten"][i] = -3e-6 * lowz * rng.uniform(0.5, 1.5)
            t3["rthdynten"][i] = 2e-4 * np.exp(-((z - 6000.0) / 2500.0) ** 2) * rng.uniform(0.5, 1.5)
            t3["rqvblten"][i] = -t3["rqvblten"][i]
        t3["rthratenlw"][i] = -2e-5 * rng.uniform(0.5, 1.5, nV)
        t3["rthratensw"][i] = 1e-5 * rng.uniform(0, 1, nV)
    coe = rng.integers(1, nCt + 1, size=(nE + 1, 2)).astype(np.int32)
    env = rng.normal(size=(nE + 1, 3))
    env /= np.linalg.norm(env, axis=1, keepdims=True)
    east = rng.normal(size=(nCt + 1, 3))
    east /= np.linalg.norm(east, axis=1, keepdims=True)
    north = rng.normal(size=(nCt + 1, 3))
    north -= (north * east).sum(1, keepdims=True) * east
    north /= np.linalg.norm(north, axis=1, keepdims=True)
    with open(path, "wb") as o:
        put(o, "hdr", np.array([nC, nCt, nEs, nE, nV, nsteps, acrain, 0, cdeep, cshal], dtype=np.int32), kind=1)
        put(o, "dtl", np.array([float(dt), float(len_disp), 2.0]))
        # mesh
        put(o, "latCell", rng.uniform(-1.5, 1.5, nCt + 1))
        put(o, "lonCell", rng.uniform(0, 6.28, nCt + 1))
        put(o, "fzm", rng.uniform(0.4, 0.6, nV))
        put(o, "fzp", rng.uniform(0.4, 0.6, nV))
        put(o, "rdzw", rng.uniform(1e-3, 1e-2, nV))
        put(o, "areaCell", rng.uniform(4e9, 6e9, nCt + 1))
        put(o, "zgrid", zg)
        put(o, "zz", zz)
        put(o, "meshDensity", rng.uniform(0.2, 1.0, nCt + 1))
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
        put(o, "surface_pressure", sp)
        for n in ["exner", "pressure_base", "pressure_p", "uReconstructZonal", "uReconstructMeridional"]:
            put(o, n, f[n])
        put(o, "rtheta_base", f["rho_zz"] * (f["theta_m"] / (1.0 + RV / RGAS * np.maximum(sc[:, :, 0], 0.0))))
        put(o, "rtheta_p", f["rho_zz"] * f["theta_m"] - f["rho_zz"] * (f["theta_m"] / (1.0 + RV / RGAS * np.maximum(sc[:, :, 0], 0.0))))
        put(o, "rho_edge", rng.uniform(0.1, 1.3, (nE + 1, nV)))
        put(o, "tend_u_phys", rng.uniform(-1, 1, (nE + 1, nV)))
        # sfc_input
        put(o, "xland", xland)
        # diag_physics
        put(o, "plrad", rng.uniform(-1, 1, nCt + 1))
        put(o, "hfx", hfx)
        put(o, "qfx", qfx)
        put(o, "kpbl", kpbl, kind=1)
        put(o, "gsw", rng.uniform(0, 900, nCt + 1))
        put(o, "cubot", rng.uniform(1, 10, nCt + 1))
        put(o, "cutop", rng.uniform(10, 40, nCt + 1))
        put(o, "xmb_total", rng.uniform(0, 1, nCt + 1))
        put(o, "xmb_shallow", rng.uniform(0, 1, nCt + 1))
        for n in ["k22_shallow", "kbcon_shallow", "ktop_shallow", "ktop_deep"]:
            put(o, n, rng.integers(0, 20, nCt + 1).astype(np.int32), kind=1)
        put(o, "qc_cu", rng.uniform(0, 1e-3, (nCt + 1, nV)))
        put(o, "qi_cu", rng.uniform(0, 1e-3, (nCt + 1, nV)))
        put(o, "cuprec", rng.uniform(0, 1e-3, nCt + 1))
        put(o, "raincv", rng.uniform(0, 1.0, nCt + 1))
        put(o, "rainc", rng.uniform(0, 4.0, nCt + 1))
        put(o, "i_rainc", rng.integers(0, 3, nCt + 1).astype(np.int32), kind=1)
        # tend_physics
        for n, a in t3.items():
            put(o, n, a)
        for n in ["rqcblten", "rqiblten", "rqsblten", "rncblten", "rniblten", "rnifablten", "rnwfablten"]:
            put(o, n, 1e-6 * rng.uniform(-1, 1, (nCt + 1, nV)))
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
    bad, total, phase = 0, 0, ""
    for (name, kind, dims, x), (_, _, _, y) in zip(a, b):
        if name.startswith("#"):
            phase = name[1:]
            continue
        total += x.size
        d = np.nonzero(x != y)[0]
        if d.size:
            bad += 1
            j = int(d[0])
            vx, vy = (x[j:j + 1].view("<f8")[0], y[j:j + 1].view("<f8")[0]) if kind == 0 else (x[j], y[j])
            idx = np.unravel_index(j, dims[::-1])[::-1]
            if bad <= 25:
                print(f"{phase} {name}: {d.size} of {x.size} differ, first at {tuple(int(t) + 1 for t in idx)}: "
                      f"ref {vx!r} cpp {vy!r}")
    print(f"{total} values in {len(a)} records, {'IDENTICAL' if not bad else str(bad) + ' records differ'}")
    return 1 if bad else 0


def stats(path, nC):
    """What the convection did in a reference output: counts over the cells 1..nC of the last step."""
    rs = records(path)
    last = {}
    for name, kind, dims, x in rs:
        if name.startswith("#"):
            continue
        last[name] = (dims, x.view("<f8") if kind == 0 else x)
    _, rc = last["raincv"]
    _, xs = last["xmb_shallow"]
    _, xt = last["xmb_total"]
    _, ks = last["k22_shallow"]
    dims2, qc = last["qc_cu"]
    _, qi = last["qi_cu"]
    qcc = qc.reshape(dims2[::-1])[:nC]
    qic = qi.reshape(dims2[::-1])[:nC]
    print(f"deep precipitation in {int((rc[:nC] > 0).sum())} of {nC} cells, xmb_total > 0 in {int((xt[:nC] > 0).sum())}; "
          f"shallow xmb > 0 in {int((xs[:nC] > 0).sum())}, k22 = kpbl branch in {int((ks[:nC] > 5).sum())}; "
          f"qc_cu > 0 in {int((qcc.max(axis=1) > 0).sum())}, qi_cu > 0 in {int((qic.max(axis=1) > 0).sum())}")


if __name__ == "__main__":
    if sys.argv[1] == "make":
        make(sys.argv[2], *map(int, sys.argv[3:7]), float(sys.argv[7]), *map(int, sys.argv[8:11]), float(sys.argv[11]))
    elif sys.argv[1] == "stats":
        stats(sys.argv[2], int(sys.argv[3]))
    else:
        sys.exit(compare(sys.argv[2], sys.argv[3]))
