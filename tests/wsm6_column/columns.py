# SPDX-License-Identifier: GPL-3.0-or-later
"""Synthetic columns for the WSM6 column check, and the comparison of the two programs' outputs.

  columns.py make <in.bin> <nCells> <nVertLevels> <nsteps> <mp_top_level> <dt> <seed>
  columns.py compare <in.bin> <ref out.bin> <cpp out.bin> [l_diags]

in.bin (native-endian stream, read by harness.F90 and test_wsm6_column.cpp): int32 nCells, nVertLevels, nsteps,
mp_top_level, nscal (6), float64 dt, then in Fortran order zgrid(nV+1,nC) zz rho_zz theta_m (nV,nC) w(nV+1,nC)
scalars(6,nV,nC) exner pressure_base pressure_p exner_base rtheta_base rtheta_p (nV,nC) surface_pressure(nC)
rainprod evapprod re_cloud re_ice re_snow (nV,nC) rainnc snownc graupelnc (nC).
out.bin: per step, the OUT fields below in that order, then (l_diags) the DIAGS fields.

Cell i gets regime i % 5: 0 warm rain (cloud and rain below the freezing level), 1 cloud ice (cold, ice
supersaturated), 2 snow, 3 graupel (deep mixed-phase column), 4 everything with some negative mixing ratios from
"the dynamics" and dry layers. Temperatures, pressures and densities are hydrostatic and consistent; amplitudes
are perturbed per cell.
"""
import sys

import numpy as np

RGAS, RV, G = 287.0, 461.6, 9.80616
CP = 7.0 * RGAS / 2.0
P0 = 1.0e5
OUT = [("scalars", "s"), ("theta_m", "v"), ("rtheta_p", "v"), ("exner", "v"), ("pressure_p", "v"),
       ("surface_pressure", "c"), ("tend_sfc_pressure", "c"), ("rt_diabatic_tend", "v"), ("dtheta_dt_mp", "v"),
       ("precipw", "c"), ("rainnc", "c"), ("rainncv", "c"), ("snownc", "c"), ("snowncv", "c"), ("graupelnc", "c"),
       ("graupelncv", "c"), ("sr", "c"), ("rainprod", "v"), ("evapprod", "v"), ("re_cloud", "v"), ("re_ice", "v"),
       ("re_snow", "v"), ("rthmpten", "v"), ("rqvmpten", "v"), ("rqcmpten", "v"), ("rqrmpten", "v"),
       ("rqimpten", "v"), ("rqsmpten", "v"), ("rqgmpten", "v")]
DIAGS = [("refl10cm", "v"), ("refl10cm_max", "c"), ("refl10cm_1km", "c"), ("refl10cm_1km_max", "c"),
         ("relhum", "v")]


def make(path, nC, nV, nsteps, mp_top, dt, seed):
    rng = np.random.default_rng(seed)
    ztop = 30000.0
    eta = np.linspace(0.0, 1.0, nV + 1)
    zg = np.empty((nC, nV + 1))
    zz = np.empty((nC, nV))
    names = ["rho_zz", "theta_m", "exner", "pressure_base", "pressure_p", "exner_base", "rtheta_base", "rtheta_p",
             "rainprod", "evapprod", "re_cloud", "re_ice", "re_snow"]
    f = {n: np.empty((nC, nV)) for n in names}
    w = np.empty((nC, nV + 1))
    sc = np.zeros((nC, nV, 6))
    sfc = np.empty(nC)
    acc = np.empty((3, nC))
    for i in range(nC):
        reg = i % 5
        hs = rng.uniform(0.0, 1500.0)
        zg[i] = hs + eta * (ztop - hs) + 30.0 * np.sin(np.pi * eta) * rng.uniform(-1, 1)
        zm = 0.5 * (zg[i, 1:] + zg[i, :-1])
        zz[i] = (ztop / (ztop - hs)) * (1.0 + 0.01 * rng.uniform(-1, 1, nV))
        tsfc = {0: 303.0, 1: 255.0, 2: 268.0, 3: 296.0, 4: 290.0}[reg] + rng.uniform(-3, 3)
        T = np.maximum(tsfc - 0.0065 * (zm - hs), 205.0 + rng.uniform(0, 5))
        p = np.empty(nV)
        ps = rng.uniform(85000.0, 102000.0) * np.exp(-hs / 8000.0)
        pk = ps
        for k in range(nV):
            dz = zm[k] - (hs if k == 0 else zm[k - 1])
            Tk = T[k] if k == 0 else 0.5 * (T[k] + T[k - 1])
            pk = pk * np.exp(-G * dz / (RGAS * Tk))
            p[k] = pk
        es = 610.78 * np.exp(17.27 * (T - 273.16) / (T - 35.86))
        qs = 0.622 * es / (p - es)
        rh = {0: 1.04, 1: 0.85, 2: 0.95, 3: 1.02, 4: 0.7}[reg] * rng.uniform(0.85, 1.08, nV)
        qv = np.clip(rh * qs, 1e-8, 0.03)
        qv[zm > 16000.0] = 3e-6
        lev = np.exp(-((zm - rng.uniform(1500, 9000)) / rng.uniform(1500, 4000)) ** 2)
        amp = rng.uniform(0.3, 1.5)
        qc = qi = qr = qsn = qg = np.zeros(nV)
        if reg in (0, 3, 4):
            qc = 1.2e-3 * amp * lev * (T > 250.0)
            qr = 2.5e-3 * amp * np.exp(-((zm - hs - 1500.0) / 2500.0) ** 2) * (T > 268.0)
        if reg in (1, 3, 4):
            qi = 2e-4 * amp * lev * (T < 268.0)
        if reg in (2, 3, 4):
            qsn = 8e-4 * amp * lev * (T < 275.0)
        if reg in (3, 4):
            qg = 2e-3 * amp * np.exp(-((zm - hs - 4000.0) / 2500.0) ** 2)
        mix = [qv, qc, qr, qi, qsn, qg]
        for n in range(6):
            sc[i, :, n] = mix[n]
        if reg == 4:  # negative values from the transport, a dry layer
            sc[i, ::7, 1:] = -1e-9 * rng.uniform(0, 1, (len(range(0, nV, 7)), 5))
            sc[i, 3, 0] = -1e-7
            sc[i, nV // 2, 1:] = 0.0
        th = T * (P0 / p) ** (RGAS / CP)
        thm = th * (1.0 + RV / RGAS * np.maximum(qv, 0.0))
        rho = p / (RGAS * T * (1.0 + 0.608 * np.maximum(qv, 0.0)))
        rho_zz = rho / zz[i]
        # base state: a slightly different column (the perturbation is not zero)
        pb = p * (1.0 + rng.uniform(-2e-3, 2e-3, nV))
        exb = (pb / P0) ** (RGAS / CP)
        rtb = rho_zz * th * (1.0 + rng.uniform(-3e-3, 3e-3, nV))
        f["rho_zz"][i] = rho_zz
        f["theta_m"][i] = thm
        f["exner"][i] = (p / P0) ** (RGAS / CP)
        f["pressure_base"][i] = pb
        f["pressure_p"][i] = p - pb
        f["exner_base"][i] = exb
        f["rtheta_base"][i] = rtb
        f["rtheta_p"][i] = rho_zz * thm - rtb
        f["rainprod"][i] = rng.uniform(0, 1e-8, nV)
        f["evapprod"][i] = rng.uniform(0, 1e-8, nV)
        f["re_cloud"][i] = 2.49e-6
        f["re_ice"][i] = 4.99e-6
        f["re_snow"][i] = 9.99e-6
        w[i] = rng.uniform(-1.0, 2.0, nV + 1)
        w[i, 0] = 0.0
        sfc[i] = ps
        acc[:, i] = rng.uniform(0.0, 5.0, 3)
    with open(path, "wb") as o:
        np.array([nC, nV, nsteps, mp_top, 6], dtype="<i4").tofile(o)
        np.array([dt], dtype="<f8").tofile(o)
        ft = lambda a: np.ascontiguousarray(a, dtype="<f8").tofile(o)  # (nC, ...) C order = Fortran (..., nC)
        ft(zg); ft(zz); ft(f["rho_zz"]); ft(f["theta_m"]); ft(w); ft(sc)
        for n in ["exner", "pressure_base", "pressure_p", "exner_base", "rtheta_base", "rtheta_p"]:
            ft(f[n])
        ft(sfc)
        for n in ["rainprod", "evapprod", "re_cloud", "re_ice", "re_snow"]:
            ft(f[n])
        ft(acc[0]); ft(acc[1]); ft(acc[2])


def compare(inp, ref, cpp, diags=None):
    out = OUT + (DIAGS if diags == "l_diags" else [])
    hdr = np.fromfile(inp, dtype="<i4", count=5)
    nC, nV, nsteps = int(hdr[0]), int(hdr[1]), int(hdr[2])
    a = np.fromfile(ref, dtype="<u8")
    b = np.fromfile(cpp, dtype="<u8")
    size = {"s": 6 * nV * nC, "v": nV * nC, "c": nC}
    per = sum(size[t] for _, t in out)
    if a.size != per * nsteps or b.size != per * nsteps:
        print(f"FAIL: sizes ref {a.size} cpp {b.size}, expected {per * nsteps}")
        return 1
    bad = 0
    off = 0
    for s in range(nsteps):
        for name, t in out:
            n = size[t]
            d = np.nonzero(a[off:off + n] != b[off:off + n])[0]
            if d.size:
                bad += 1
                j = d[0]
                x, y = a[off + j:off + j + 1].view("<f8")[0], b[off + j:off + j + 1].view("<f8")[0]
                print(f"step {s + 1} {name}: {d.size} of {n} differ, first at flat index {j}: ref {x!r} cpp {y!r}")
            off += n
    hc = np.frombuffer(a.tobytes(), dtype="<f8")
    print(f"{nC} cells x {nV} levels, {nsteps} steps: {per * nsteps} values, "
          f"{'IDENTICAL' if not bad else str(bad) + ' field-steps differ'}; nonzero values {np.count_nonzero(hc)}")
    return 1 if bad else 0


if __name__ == "__main__":
    if sys.argv[1] == "make":
        make(sys.argv[2], *map(int, sys.argv[3:7]), float(sys.argv[7]), int(sys.argv[8]))
    else:
        sys.exit(compare(*sys.argv[2:6]))
