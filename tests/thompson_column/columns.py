# SPDX-License-Identifier: GPL-3.0-or-later
"""Synthetic columns for the Thompson column check, and the comparison of the two programs' outputs.

  columns.py make <in.bin> <nCells> <nVertLevels> <nsteps> <mp_top_level> <dt> <seed>
  columns.py compare <in.bin> <ref out.bin> <cpp out.bin> [l_diags]
  columns.py compare_init <ref init.bin> <cpp init.bin>

in.bin (little-endian stream, read by harness.F90 and test_thompson_column.cpp): int32 nCells, nVertLevels, nsteps,
mp_top_level, nscal (8: qv qc qr qi qs qg ni nr), float64 dt, then in Fortran order zgrid(nV+1,nC) zz rho_zz theta_m
(nV,nC) w(nV+1,nC) scalars(8,nV,nC) exner pressure_base pressure_p exner_base rtheta_base rtheta_p (nV,nC)
surface_pressure(nC) rainprod evapprod re_cloud re_ice re_snow (nV,nC) rainnc snownc graupelnc landmask nt_c mu_c
(nC; landmask as a double 0, 1 or 2).
out.bin: nt_c, mu_c (after init_thompson_clouddroplets_forMPAS), dt_microp, n_microp, then per step the OUT fields
below in that order, then (l_diags) the DIAGS fields.
init.bin: thompson_init's module variables (dump.c / dump_thompson_init), compared record by record.

Cell i gets regime i % 7: 0 warm rain (cloud and rain below the freezing level), 1 cloud ice (cold, ice
supersaturated), 2 snow, 3 graupel (deep mixed-phase column), 4 everything with some negative mixing ratios from
"the dynamics" and dry layers, 5 clear and subsaturated (mp_thompson returns at once), 6 clear and supersaturated
near the ground (condensation starts). Temperatures, pressures and densities are hydrostatic and consistent;
amplitudes are perturbed per cell. ni and nr are positive where qi / qr are, with some zero or tiny numbers.
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
       ("rqimpten", "v"), ("rqsmpten", "v"), ("rqgmpten", "v"), ("rnimpten", "v"), ("rnrmpten", "v")]
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
    sc = np.zeros((nC, nV, 8))
    lmask = np.empty(nC)
    ntc0 = np.empty(nC)
    muc0 = np.empty(nC)
    sfc = np.empty(nC)
    acc = np.empty((3, nC))
    for i in range(nC):
        reg = i % 7
        hs = rng.uniform(0.0, 1500.0)
        zg[i] = hs + eta * (ztop - hs) + 30.0 * np.sin(np.pi * eta) * rng.uniform(-1, 1)
        zm = 0.5 * (zg[i, 1:] + zg[i, :-1])
        zz[i] = (ztop / (ztop - hs)) * (1.0 + 0.01 * rng.uniform(-1, 1, nV))
        tsfc = {0: 303.0, 1: 255.0, 2: 268.0, 3: 296.0, 4: 290.0, 5: 285.0, 6: 280.0}[reg] + rng.uniform(-3, 3)
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
        rh = {0: 1.04, 1: 0.85, 2: 0.95, 3: 1.02, 4: 0.7, 5: 0.5, 6: 0.9}[reg] * rng.uniform(0.85, 1.08, nV)
        if reg == 6:
            rh[: max(2, nV // 6)] = 1.03
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
        ni = np.where(qi > 1e-12, 10.0 ** rng.uniform(3.0, 5.5, nV), 0.0)
        nr = np.where(qr > 1e-12, 10.0 ** rng.uniform(2.5, 5.0, nV), 0.0)
        ni[::11] = np.where(qi[::11] > 0, 1e-9, 0.0)  # tiny numbers: the R2 limits
        nr[::13] = np.where(qr[::13] > 0, 0.0, 0.0)
        mix = [qv, qc, qr, qi, qsn, qg, ni, nr]
        for n in range(8):
            sc[i, :, n] = mix[n]
        lmask[i] = i % 3
        ntc0[i] = 150.e6
        muc0[i] = 9.0
        if reg == 4:  # negative values from the transport, a dry layer
            sc[i, ::7, 1:6] = -1e-9 * rng.uniform(0, 1, (len(range(0, nV, 7)), 5))
            sc[i, 3, 0] = -1e-7
            sc[i, nV // 2, 1:] = 0.0
            sc[i, 5::9, 6:] = -1.0
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
        np.array([nC, nV, nsteps, mp_top, 8], dtype="<i4").tofile(o)
        np.array([dt], dtype="<f8").tofile(o)
        ft = lambda a: np.ascontiguousarray(a, dtype="<f8").tofile(o)  # (nC, ...) C order = Fortran (..., nC)
        ft(zg); ft(zz); ft(f["rho_zz"]); ft(f["theta_m"]); ft(w); ft(sc)
        for n in ["exner", "pressure_base", "pressure_p", "exner_base", "rtheta_base", "rtheta_p"]:
            ft(f[n])
        ft(sfc)
        for n in ["rainprod", "evapprod", "re_cloud", "re_ice", "re_snow"]:
            ft(f[n])
        ft(acc[0]); ft(acc[1]); ft(acc[2]); ft(lmask); ft(ntc0); ft(muc0)


def compare(inp, ref, cpp, diags=None):
    out = OUT + (DIAGS if diags == "l_diags" else [])
    hdr = np.fromfile(inp, dtype="<i4", count=5)
    nC, nV, nsteps = int(hdr[0]), int(hdr[1]), int(hdr[2])
    a = np.fromfile(ref, dtype="<u8")
    b = np.fromfile(cpp, dtype="<u8")
    size = {"s": 8 * nV * nC, "v": nV * nC, "c": nC}
    per = sum(size[t] for _, t in out)
    head = 2 * nC + 2
    if a.size != head + per * nsteps or b.size != head + per * nsteps:
        print(f"FAIL: sizes ref {a.size} cpp {b.size}, expected {head + per * nsteps}")
        return 1
    bad = 0
    for name, o, n in (("nt_c", 0, nC), ("mu_c", nC, nC), ("dt_microp/n_microp", 2 * nC, 2)):
        d = np.nonzero(a[o:o + n] != b[o:o + n])[0]
        if d.size:
            bad += 1
            print(f"{name}: {d.size} of {n} differ")
    off = head
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
    print(f"{nC} cells x {nV} levels, {nsteps} steps: {head + per * nsteps} values, "
          f"{'IDENTICAL' if not bad else str(bad) + ' field-steps differ'}; nonzero values {np.count_nonzero(hc)}")
    return 1 if bad else 0


def read_dump(path):
    raw = open(path, "rb").read()
    out, o = [], 0
    while o < len(raw):
        name = raw[o:o + 32].rstrip(b"\0").decode()
        n = int(np.frombuffer(raw, dtype="<i8", count=1, offset=o + 32)[0])
        o += 40
        cnt = 2 if n == -1 else n
        out.append((name, n, raw[o:o + 8 * cnt]))
        o += 8 * cnt
    return out


def compare_init(ref, cpp):
    a, b = read_dump(ref), read_dump(cpp)
    if [x[:2] for x in a] != [x[:2] for x in b]:
        print("FAIL: the two dumps have different records")
        return 1
    bad, nval = 0, 0
    for (name, n, x), (_, _, y) in zip(a, b):
        nval += abs(n) if n != -1 else int(np.frombuffer(x, "<f8")[0])
        if x != y:
            bad += 1
            if n == -1:
                print(f"{name}: table hash differs")
            else:
                u, v = np.frombuffer(x, "<u8"), np.frombuffer(y, "<u8")
                d = np.nonzero(u != v)[0]
                print(f"{name}: {d.size} of {n} differ, first at {d[0]}: ref {np.frombuffer(x, '<f8')[d[0]]!r} "
                      f"cpp {np.frombuffer(y, '<f8')[d[0]]!r}")
    print(f"thompson_init: {len(a)} records, {nval} values: {'IDENTICAL' if not bad else str(bad) + ' records differ'}")
    return 1 if bad else 0


if __name__ == "__main__":
    if sys.argv[1] == "compare_init":
        sys.exit(compare_init(sys.argv[2], sys.argv[3]))
    if sys.argv[1] == "make":
        make(sys.argv[2], *map(int, sys.argv[3:7]), float(sys.argv[7]), int(sys.argv[8]))
    else:
        sys.exit(compare(*sys.argv[2:6]))
