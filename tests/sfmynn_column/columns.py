# SPDX-License-Identifier: GPL-3.0-or-later
"""Synthetic columns for the MYNN surface layer column check, and the comparison of the two programs' outputs.

  columns.py make <in.bin> <nCellsSolve> <nVertLevels> <nsteps> <seed> <frac_seaice> <do_restart> <len_disp>
  columns.py compare <in.bin> <ref out.bin> <cpp out.bin>   (the input gives nsteps: every step must be in both)

Both files are native-endian streams of records: name (32 bytes, blank padded), int32 kind (0 float64, 1 int32),
int32 ndim, int32 dims[4], then the data in Fortran order. in.bin starts with the records "hdr" (int32: nCellsSolve,
nCells, nVertLevels, nsteps, frac_seaice, do_restart, 0, 0) and "dtl" (float64: dt, config_len_disp), then the pool
arrays (their Registry names; the harnesses register each in its pool). The outputs are, in order: the four psi
tables of sf_mynn_init; then per step the snapshots "A" (after MPAS_to_physics: every pool array and the
mpas_atmphys_vars arrays the surface layer reads, names starting with "@") and "B" (after driver_sfclayer: every pool
array, then, with config_frac_seaice, the *_sea arrays that outlive the driver), each a marker record first.

nCells = nCellsSolve + 2 (two halo cells MPAS_to_physics and the driver never touch). Cell i (1-based) gets regime
(i - 1) % 10:
  0 land, daytime convective (hot ground, light wind; qsfc = 0 for every third cell: the water formula); by
    (i - 1) // 10 % 3: pblh near 3 km (wstar's 4000 m cap), snowh 0.09 (just below the snow threshold), or a thin
    (za = 3 m), very rough (z0 ~ 2 m), calm, strongly heated layer (the 2 m theta fallback)
  1 land, stable night, calm (wind down to 0.02 m/s: wspd = wmin with len_disp 3 km; ust = max(ust, 0.005))
  2 sea, unstable (warm sea, moist): charnock_1955 and fairall_etal_2003 (restar <= 2 and > 2)
  3 sea, stable, strong wind (wsp10m > 18: the Charnock coefficient at its cap)
  4 sea ice: xice in [0.3, 1], xland 1 (as MPAS sets it) or 2, ice temperatures (tsk < 273.15: the ice formula),
    snow on the ice for half of them (andreas_2002 in the first call, the open-water call over the rest)
  5 land, snow, strong wind and rough (andreas_2002 with ren2 > 1000 once ust has grown)
  6 land, snow, calm and smooth (andreas_2002 with 0.135 < ren2 < 2.5)
  7 sea, cold air over warm sea, calm, shallow pblh (br at its lower limits, z/L beyond -10: the full unstable
    functions)
  8 land, strong surface inversion, calm (z/L at its cap 20: the full stable functions; every other such cell has
    za = 20.01 + 2 z0, so that (10 + z0)/L falls in 9.99..10, the end of MPAS's tables: nzol + 1 = 1000)
  9 sea, xice 0.01 (below the threshold) or 0.02 (on it), cloud water at the lowest level, ust < 0.001
The first two layers vary with the cell ((7 i + seed) % 4): za = 20 m (the flux-profile 10 m wind), 10 m (the
neutral log), 5 m with za2 in (7, 13) (the second level's wind) and 3 m with za2 >= 13 (the neutral log again).
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
    """The records of a column file; a short header or short data (a truncated file) is an error."""
    out = []
    with open(path, "rb") as f:
        while True:
            h = f.read(NAMELEN)
            if not h:
                break
            raw = f.read(24)
            if len(h) != NAMELEN or len(raw) != 24:
                raise ValueError(f"{path}: record {len(out) + 1} has a short header (truncated file)")
            name = h.decode().strip()
            meta = np.frombuffer(raw, dtype="<i4")
            kind, nd, dims = int(meta[0]), int(meta[1]), [int(x) for x in meta[2:6]]
            if kind not in (0, 1) or not 1 <= nd <= 4 or min(dims[:nd]) < 0:
                raise ValueError(f"{path}: record {len(out) + 1} ({name!r}) has a bad header {meta.tolist()}")
            n = int(np.prod(dims[:nd]))
            size = n * (8 if kind == 0 else 4)
            body = f.read(size)
            if len(body) != size:
                raise ValueError(f"{path}: record {len(out) + 1} ({name}) has {len(body)} of {size} bytes "
                                 f"(truncated file)")
            data = np.frombuffer(body, dtype="<u8" if kind == 0 else "<i4")
            out.append((name, kind, dims[:nd], data))
    return out


def layout_problems(inp, recs):
    """What a complete output of the input's nsteps steps holds: every step's phase markers (A, B), each followed
    by data records, and at least one value."""
    hdr = [x for n, _, _, x in records(inp) if n == "hdr"]
    if len(hdr) != 1:
        return [f"{inp}: no hdr record"]
    nsteps = int(hdr[0][3])
    if nsteps < 1:
        return [f"{inp}: nsteps = {nsteps}"]
    want = {(ph, s) for ph in ('A', 'B') for s in range(1, nsteps + 1)}
    seen, problems = set(), []
    for i, (name, kind, dims, x) in enumerate(recs):
        if name.startswith("#"):
            if kind != 1 or x.size != 1:
                problems.append(f"marker {name} is not one integer")
                continue
            seen.add((name[1:], int(x[0])))
            if i + 1 == len(recs) or recs[i + 1][0].startswith("#"):
                problems.append(f"marker {name} {int(x[0])} is followed by no data record")
    unknown = sorted(seen - want)
    if unknown:
        problems.append(f"{len(unknown)} phase markers not of steps 1..{nsteps} (first {unknown[0]})")
    missing = sorted(want - seen)
    if missing:
        problems.append(f"{len(missing)} phase markers missing (first {missing[0]}): the run did not finish")
    if not any(x.size for n, _, _, x in recs if not n.startswith("#")):
        problems.append("no value written")
    return problems


LAND = (0, 1, 5, 6, 8)


def make(path, nC, nV, nsteps, seed, seaice, restart, len_disp):
    rng = np.random.default_rng(seed)
    nCt = nC + 2
    ztop = 25000.0
    shape_c = (nCt + 1,)
    zg = np.empty((nCt + 1, nV + 1))
    zz = np.empty((nCt + 1, nV))
    names3 = ["rho_zz", "theta_m", "exner", "pressure_base", "pressure_p", "rtheta_base", "rtheta_p",
              "uReconstructZonal", "uReconstructMeridional"]
    f = {n: np.empty((nCt + 1, nV)) for n in names3}
    w = np.empty((nCt + 1, nV + 1))
    sc = np.zeros((nCt + 1, nV, 6))
    sfc = {n: np.empty(shape_c) for n in ["skintemp", "xland", "sst", "xice", "snowh"]}
    d1 = {n: np.empty(shape_c) for n in
          ["plrad", "hpbl", "mavail", "br", "cpm", "chs", "chs2", "cqs2", "flhc", "flqc", "gz1oz0", "hfx", "qfx",
           "qgh", "qsfc", "lh", "mol", "psih", "psim", "regime", "rmol", "ust", "ustm", "wspd", "znt", "zol", "q2",
           "t2m", "th2m", "u10", "v10", "cd", "cda", "ck", "cka", "fh", "fm", "ch", "qcg"]}
    surface_pressure = np.empty(shape_c)
    for i in range(nCt + 1):
        reg = i % 10  # array index i = cell i + 1
        land = reg in LAND
        znt8 = rng.uniform(0.01, 0.1) if reg == 8 else 0.0
        hs = rng.uniform(0.0, 1200.0) if land else 0.0
        dz1, dz2 = [(40.0, 45.0), (20.0, 25.0), (10.0, 5.0), (6.0, 30.0)][(7 * i + seed) % 4]
        dz1 *= rng.uniform(0.98, 1.02)
        thin = reg == 0 and (i // 10) % 3 == 2  # thin first layer, very rough, calm, strongly heated
        if thin:
            dz1, dz2 = 6.0, 30.0
        if reg == 8 and (i // 10) % 2 == 0:  # za = 20.01 + 2 z0: with z/L at its cap 20, (10 + z0)/L is 9.995
            dz1, dz2 = 40.02 + 4.0 * znt8, 45.0
        z2 = dz1 + dz2
        eta = np.linspace(0.0, 1.0, nV - 1) ** 1.5
        zrel = np.concatenate(([0.0, dz1], z2 + eta * (ztop - hs - z2)))
        zg[i] = hs + zrel
        zm = 0.5 * (zg[i, 1:] + zg[i, :-1])
        zz[i] = (ztop / (ztop - hs)) * (1.0 + 0.003 * rng.uniform(-1, 1, nV))
        z = zm - hs
        tsfc = {0: 302.0, 1: 281.0, 2: 299.0, 3: 284.0, 4: 258.0, 5: 266.0, 6: 262.0, 7: 279.0, 8: 274.0,
                9: 289.0}[reg] + rng.uniform(-2, 2)
        lapse = np.full(nV, 0.0065)
        if reg in (0, 2, 7):  # well mixed (dry adiabatic) below the PBL top
            lapse[z < rng.uniform(800, 1800)] = 0.0098
        if reg in (1, 3, 4, 6, 8):  # surface inversion
            lapse[z < rng.uniform(150, 400)] = -0.012 * rng.uniform(0.5, 1.5) * (2.5 if reg == 8 else 1.0)
        if reg == 9:  # mixed layer capped by an inversion at ~1 km
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
        rh = {0: 0.45, 1: 0.8, 2: 0.85, 3: 0.9, 4: 0.85, 5: 0.7, 6: 0.9, 7: 0.6, 8: 0.75, 9: 0.97}[reg]
        qv = np.clip(rh * rng.uniform(0.9, 1.05, nV) * qs, 1e-8, 0.03)
        qv[zm > 12000.0] = 3e-6
        qc = np.zeros(nV)
        if reg == 9:
            qc[(z > zi - 300) & (z < zi)] = rng.uniform(2e-4, 6e-4)
            qc[0] = 1e-5 * rng.uniform(0, 1)
        qr = 1e-5 * rng.uniform(0, 1, nV) * (z < 3000)
        qsn = 1e-5 * rng.uniform(0, 1, nV) * (T < 270)
        mix = [qv, qc, qr, np.zeros(nV), qsn, np.zeros(nV)]
        for n in range(6):
            sc[i, :, n] = mix[n]
        th = T * (P0 / p) ** (RGAS / CP)
        thm = th * (1.0 + RV / RGAS * qv)
        rho = p / (RGAS * T * (1.0 + RV / RGAS * qv))  # dry density
        pb = p * (1.0 + rng.uniform(-2e-3, 2e-3, nV))
        f["rho_zz"][i] = rho / zz[i]
        f["theta_m"][i] = thm
        f["exner"][i] = (p / P0) ** (RGAS / CP)
        f["pressure_base"][i] = pb
        f["pressure_p"][i] = p - pb
        f["rtheta_base"][i] = f["rho_zz"][i] * th
        f["rtheta_p"][i] = f["rho_zz"][i] * thm - f["rtheta_base"][i]
        spd = {0: rng.uniform(1.0, 4.0), 1: rng.uniform(0.02, 1.5), 2: rng.uniform(3.0, 9.0),
               3: rng.uniform(18.0, 26.0), 4: rng.uniform(4.0, 12.0), 5: rng.uniform(18.0, 28.0),
               6: rng.uniform(0.3, 1.0), 7: rng.uniform(0.2, 0.8), 8: rng.uniform(0.3, 1.5), 9: rng.uniform(4.0, 9.0)}[reg]
        if thin:
            spd *= 0.05
        prof = np.log(1.0 + z / 0.1) / np.log(1.0 + 1000.0 / 0.1)
        dirn = rng.uniform(0, 2 * np.pi) + 0.3 * z / 5000.0
        f["uReconstructZonal"][i] = spd * np.minimum(prof, 1.5 + z / 10000.0) * np.cos(dirn)
        f["uReconstructMeridional"][i] = spd * np.minimum(prof, 1.5 + z / 10000.0) * np.sin(dirn)
        w[i] = 0.05 * rng.uniform(-1, 1, nV + 1)
        w[i, 0] = 0.0
        surface_pressure[i] = ps * rng.uniform(0.999, 1.001)
        # surface: the ground, sea or ice temperature relative to the air
        dT = {0: 8.0, 1: -6.0, 2: 1.5, 3: -2.5, 4: -3.0, 5: -2.0, 6: -4.0, 7: 15.0, 8: -14.0, 9: 0.8}[reg]
        if thin:
            dT = 25.0
        sfc["skintemp"][i] = T[0] + 0.0098 * z[0] + dT * rng.uniform(0.6, 1.4)
        sfc["xland"][i] = 1.0 if land else 2.0
        sfc["sst"][i] = sfc["skintemp"][i] if not land else 285.0
        sfc["xice"][i] = 0.0
        sfc["snowh"][i] = 0.0
        if reg == 4:
            sfc["sst"][i] = 271.0 + rng.uniform(-1, 1.5)
            sfc["xice"][i] = rng.uniform(0.3, 1.0)
            sfc["xland"][i] = 1.0 if (i // 10) % 3 != 2 else 2.0
            sfc["snowh"][i] = 0.2 if i % 2 else 0.05
        if reg == 9:
            sfc["xice"][i] = [0.01, 0.02][i % 2]
        if reg == 5:
            sfc["snowh"][i] = rng.uniform(0.2, 1.0)
        if reg == 6:
            sfc["snowh"][i] = rng.uniform(0.1, 0.3)
        if reg == 0 and (i // 10) % 3 == 1:
            sfc["snowh"][i] = 0.09  # just below the snow threshold
        # the diag_physics state the previous physics step left
        uns = reg in (0, 2, 7, 9)
        d1["hpbl"][i] = rng.uniform(800, 1800) if uns else rng.uniform(50, 300)
        if reg == 0 and (i // 10) % 3 == 0:
            d1["hpbl"][i] = rng.uniform(2800, 3500)  # min(1.5 pblh, 4000.) at its cap
        if reg == 7 or thin:
            d1["hpbl"][i] = rng.uniform(30, 80)  # small wstar: br below -4 after the first step (len_disp 3 km)
        d1["mavail"][i] = 1.0 if not land else rng.uniform(0.1, 0.5)
        d1["br"][i] = -0.5 * rng.uniform(0, 1) if uns else 0.2 * rng.uniform(0, 1)
        d1["cpm"][i] = CP * (1 + 0.8 * qv[0])
        d1["znt"][i] = {5: rng.uniform(0.5, 1.0), 6: rng.uniform(0.001, 0.005), 8: znt8}.get(
            reg, rng.uniform(0.05, 0.8) if land else rng.uniform(1e-4, 5e-4))
        d1["ust"][i] = {5: rng.uniform(1.2, 2.0), 6: rng.uniform(0.03, 0.08), 9: 0.0005}.get(reg, rng.uniform(0.1, 0.6))
        if thin:
            d1["znt"][i] = 2.0 + 0.5 * np.sin(i)
            d1["ust"][i] = 0.01
        d1["ustm"][i] = d1["ust"][i] * rng.uniform(0.9, 1.1)
        d1["mol"][i] = -rng.uniform(0.05, 0.5) if uns else rng.uniform(0.01, 0.2)
        d1["hfx"][i] = rng.uniform(50, 350) if uns else -rng.uniform(5, 60)
        d1["qfx"][i] = rng.uniform(2e-5, 2e-4) if uns else rng.uniform(-1e-6, 2e-5)
        d1["qsfc"][i] = qv[0] * rng.uniform(1.0, 1.3) if (land and i % 3) else 0.0
        d1["regime"][i] = 4.0 if uns else 1.0
        d1["wspd"][i] = spd * 0.4
        for n in ["chs", "chs2", "cqs2", "flhc", "flqc", "gz1oz0", "qgh", "lh", "psih", "psim", "rmol", "zol", "q2",
                  "t2m", "th2m", "cd", "cda", "ck", "cka", "plrad", "fh", "fm", "ch"]:
            d1[n][i] = rng.uniform(-1, 1)
        d1["qcg"][i] = 1e-4 * rng.uniform(0, 1)
        d1["u10"][i] = f["uReconstructZonal"][i, 0] * 0.8
        d1["v10"][i] = f["uReconstructMeridional"][i, 0] * 0.8
    with open(path, "wb") as o:
        put(o, "hdr", np.array([nC, nCt, nV, nsteps, seaice, restart, 0, 0], dtype=np.int32), kind=1)
        put(o, "dtl", np.array([600.0, float(len_disp)]))
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
        # sfc_input
        for n in ["skintemp", "xland", "sst", "xice", "snowh"]:
            put(o, n, sfc[n])
        # diag_physics
        for n, a in d1.items():
            put(o, n, a)


def compare(inp, ref, cpp):
    try:
        a, b = records(ref), records(cpp)
    except ValueError as exc:
        print(f"FAIL: {exc}")
        return 1
    for label, recs in (("ref", a), ("cpp", b)):
        problems = layout_problems(inp, recs)
        for x in problems:
            print(f"FAIL: {label}: {x}")
        if problems:
            return 1
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
            phase = name[1:] + str(int(x[0]))
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
        make(sys.argv[2], *map(int, sys.argv[3:10]))
    else:
        sys.exit(compare(sys.argv[2], sys.argv[3], sys.argv[4]))
