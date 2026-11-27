# SPDX-License-Identifier: GPL-3.0-or-later
"""Synthetic columns for the radiation column check, and the comparison of the two programs' outputs.

  columns.py make <in.bin> <nCellsSolve> <nVertLevels> <nsteps> <seed> <cldovrlp> <dcorrlen> <o3climatology>
  columns.py compare <in.bin> <ref out.bin> <cpp out.bin>   (the input gives nsteps: every step must be in both)

cldovrlp: config_radt_cld_overlap 0 'none', 1 'random', 2 'maximum_random' (mesoscale_reference), 3 'maximum',
5 'exponential_random'; dcorrlen: config_radt_cld_dcorrlen 0 'constant', 1 'latitude_varying'; o3climatology:
config_o3climatology 1 (the default) or 0.

Both files are native-endian streams of records: name (32 bytes, blank padded), int32 kind (0 float64, 1 int32),
int32 ndim, int32 dims[4], then the data in Fortran order. in.bin starts with "hdr" (int32: nCellsSolve,
nVertLevels, nsteps, nOznLevels, nMonths, cldovrlp, dcorrlen, o3climatology, 0, 0), "dtl" (float64: dt_radtsw,
dt_radtlw in s, config_len_disp) and "time" (float64 (4, nsteps): xtime_s, julday, curr_julday, gmt per step), then
the pool arrays (their Registry names) and the mpas_atmphys_vars arrays MPAS_to_physics fills (names starting with
"@", (nCellsSolve, levels)). The outputs are a snapshot "I" after the initialisation (pin, ozmixm) and one per step
"S" after the cloudiness, SW and LW drivers: a marker record, then every pool array.

nCells = nCellsSolve + 1 (one halo cell: init_o3climatology fills it), arrays with nCells + 1 elements as MPAS
allocates them. Cell i (1-based) gets regime (i - 1) % 8 (the cloud water, ice and snow profiles; the cloud fraction
of calc_cldfraction follows from them and the humidity):
  0 clear, dry                         1 low liquid cloud (stratocumulus), land
  2 deep mixed-phase cloud, sea        3 cirrus (ice and snow), supersaturated (rhum >= 1: cloud fraction 1)
  4 multi-layer cloud with clear gaps  5 negative and -0.0 mixing ratios (transport), tiny qv
  6 sea ice / snow-covered land, polar 7 fog at the ground, high cloud tops
Latitudes from pole to pole and longitudes around the globe give day and night columns and low suns at every step;
the model top lies at 18-32 km, so the top pressure (3-80 hPa) gives 1-20 LW buffer layers; temperatures,
pressures and densities are hydrostatic and consistent.
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
    """What a complete output of the input's nsteps steps holds: every step's phase markers, each followed by data
    records, and at least one value."""
    hdr = [x for n, _, _, x in records(inp) if n == "hdr"]
    if len(hdr) != 1:
        return [f"{inp}: no hdr record"]
    nsteps = int(hdr[0][2])
    if nsteps < 1:
        return [f"{inp}: nsteps = {nsteps}"]
    want = {(ph, s) for ph in ("S",) for s in range(1, nsteps + 1)}
    want.add(("I", 0))  # the initialisation
    seen, problems = set(), []
    for i, (name, kind, dims, x) in enumerate(recs):
        if name.startswith("#"):
            if kind != 1 or x.size != 1:
                problems.append(f"marker {name} is not one integer")
                continue
            seen.add((name[1:], int(x[0])))
            if i + 1 == len(recs) or recs[i + 1][0].startswith("#"):
                problems.append(f"marker {name} {int(x[0])} is followed by no data record")
    missing = sorted(want - seen)
    if missing:
        problems.append(f"{len(missing)} phase markers missing (first {missing[0]}): the run did not finish")
    if not any(x.size for n, _, _, x in recs if not n.startswith("#")):
        problems.append("no value written")
    return problems


def make(path, nC, nV, nsteps, seed, cldovrlp, dcorrlen, o3clim):
    rng = np.random.default_rng(seed)
    nO, nM = 59, 12
    ncell = nC + 2  # nCells + 1 elements
    lat = np.empty(ncell)
    lon = np.empty(ncell)
    sfc = {n: np.empty(ncell) for n in ["skintemp", "snow", "xice", "xland", "sfc_albedo", "sfc_emiss"]}
    pres = np.zeros((nC, nV)); pres2 = np.zeros((nC, nV + 1)); pi = np.zeros((nC, nV)); T = np.zeros((nC, nV))
    T2 = np.zeros((nC, nV + 1)); dz = np.zeros((nC, nV)); qv = np.zeros((nC, nV)); qc = np.zeros((nC, nV))
    qi = np.zeros((nC, nV)); qs = np.zeros((nC, nV))
    for i in range(ncell):
        reg = i % 8
        if reg == 6:
            lat[i] = np.deg2rad(rng.choice([-1, 1]) * rng.uniform(62, 89.5))
        else:
            lat[i] = np.deg2rad(rng.uniform(-80, 80))
        lon[i] = rng.uniform(0, 2 * np.pi)
        land = reg in (0, 1, 5, 6) and not (reg == 6 and i % 16 == 6)
        sfc["xland"][i] = 1.0 if land else 2.0
        sfc["xice"][i] = rng.uniform(0.3, 1.0) if (reg == 6 and not land) else (0.0 if rng.uniform() < 0.8 else 0.01)
        sfc["snow"][i] = rng.uniform(5.0, 300.0) if (reg == 6 and land) else 0.0
        sfc["skintemp"][i] = {6: 250.0}.get(reg, 285.0) + rng.uniform(-15, 20)
        sfc["sfc_albedo"][i] = rng.uniform(0.5, 0.85) if reg == 6 else rng.uniform(0.06, 0.3)
        sfc["sfc_emiss"][i] = rng.uniform(0.9, 1.0)
        if i >= nC:
            continue
        hs = rng.uniform(0.0, 2500.0) if land else 0.0
        ztop = rng.uniform(18000.0, 32000.0)
        eta = np.linspace(0.0, 1.0, nV + 1) ** 1.5
        zw = hs + eta * (ztop - hs)
        zm = 0.5 * (zw[1:] + zw[:-1])
        dz[i] = zw[1:] - zw[:-1]
        tsfc = sfc["skintemp"][i] + rng.uniform(-3, 3)
        tropo = rng.uniform(9000.0, 17000.0)
        tk = np.where(zm < tropo, tsfc - 0.0065 * (zm - hs), tsfc - 0.0065 * (tropo - hs) + 0.0015 * (zm - tropo))
        tk = np.maximum(tk, 190.0 + rng.uniform(0, 10))
        T[i] = tk
        ps = rng.uniform(98000.0, 103500.0) * np.exp(-hs / 8000.0)
        pw = np.empty(nV + 1)
        pw[0] = ps
        for k in range(nV):
            pw[k + 1] = pw[k] * np.exp(-G * dz[i, k] / (RGAS * tk[k]))
        pres2[i] = pw
        pres[i] = np.sqrt(pw[1:] * pw[:-1])
        pi[i] = (pres[i] / P0) ** (RGAS / CP)
        T2[i, 1:nV] = 0.5 * (tk[1:] + tk[:-1])
        T2[i, 0] = tsfc + rng.uniform(-1, 1)
        T2[i, nV] = tk[-1] + 0.5 * (tk[-1] - tk[-2])
        es = 610.78 * np.exp(17.27 * (tk - 273.16) / (tk - 35.86))
        qsat = 0.622 * es / (pres[i] - es)
        rh = np.clip(rng.uniform(0.5, 0.9) * np.exp(-(zm - hs) / 4000.0), 0.02, 1.0)
        qv[i] = rh * qsat
        z = zm - hs
        if reg == 1:  # stratocumulus
            zb = rng.uniform(400, 1200)
            m = (z > zb) & (z < zb + 600)
            qc[i, m] = rng.uniform(1e-5, 6e-4, m.sum())
            qv[i, m] = qsat[m] * rng.uniform(0.95, 1.0, m.sum())
        elif reg == 2:  # deep mixed phase
            m = (z > 800) & (z < tropo)
            qc[i, m & (tk > 253)] = rng.uniform(1e-6, 1e-3, (m & (tk > 253)).sum())
            qi[i, m & (tk < 268)] = rng.uniform(1e-7, 2e-4, (m & (tk < 268)).sum())
            qs[i, m & (tk < 273)] = rng.uniform(1e-7, 5e-4, (m & (tk < 273)).sum())
            qv[i, m] = qsat[m] * rng.uniform(0.85, 1.02, m.sum())
        elif reg == 3:  # cirrus, supersaturated
            m = (z > tropo - 3000) & (z < tropo - 300)
            qi[i, m] = rng.uniform(1e-7, 1e-4, m.sum())
            qs[i, m] = rng.uniform(0, 5e-5, m.sum())
            qv[i, m] = qsat[m] * rng.uniform(0.98, 1.1, m.sum())
        elif reg == 4:  # layers with gaps
            for zb, zt in ((500, 1500), (3000, 4500), (7000, 9000)):
                m = (z > zb) & (z < zt) & (rng.uniform(size=nV) < 0.7)
                qc[i, m & (tk > 258)] = rng.uniform(1e-6, 4e-4, (m & (tk > 258)).sum())
                qi[i, m & (tk < 263)] = rng.uniform(1e-7, 1e-4, (m & (tk < 263)).sum())
                qv[i, m] = qsat[m] * rng.uniform(0.7, 1.0, m.sum())
        elif reg == 5:  # transport noise
            qc[i] = rng.uniform(-1e-6, 2e-5, nV)
            qi[i] = rng.uniform(-1e-6, 1e-5, nV)
            qs[i] = rng.uniform(-1e-6, 1e-5, nV)
            qc[i, ::7] = -0.0
            qv[i, -5:] = rng.uniform(-1e-9, 1e-13, 5)
        elif reg == 6:  # polar low cloud
            m = (z > 200) & (z < 2000)
            qc[i, m] = rng.uniform(0, 1e-4, m.sum())
            qi[i, m] = rng.uniform(0, 5e-5, m.sum())
        elif reg == 7:  # fog and high tops
            qc[i, :3] = rng.uniform(1e-4, 5e-4, 3)
            qv[i, :3] = qsat[:3] * 1.001
            m = (z > tropo - 1500) & (z < tropo + 500)
            qi[i, m] = rng.uniform(1e-6, 5e-5, m.sum())
    # times: start 2026-10-12_00 (julday 285, gmt 0); steps 0, 6 h, 12 h, 18 h, 25 h, ... (day of year crossing)
    tm = np.zeros((nsteps, 4))
    for s in range(nsteps):
        xt = [0.0, 21600.0, 43200.0, 64800.0, 90000.0, 1800.0][s % 6] + 86400.0 * (s // 6) * 3
        jd = 285 + int(xt // 86400.0)
        utc_h = (xt % 86400.0) / 3600.0
        tm[s] = [xt, jd, (jd - 1) + utc_h / 24.0, 0.0]
    with open(path, "wb") as o:
        put(o, "hdr", np.array([nC, nV, nsteps, nO, nM, cldovrlp, dcorrlen, o3clim, 0, 0]), kind=1)
        put(o, "dtl", np.array([1800.0, 1800.0, 240000.0]))
        put(o, "time", tm)
        put(o, "latCell", lat); put(o, "lonCell", lon)
        put(o, "areaCell", np.full(ncell, 1.0e10)); put(o, "meshDensity", np.ones(ncell))
        for n in ["skintemp", "snow", "xice", "xland", "sfc_albedo", "sfc_emiss"]:
            put(o, n, sfc[n])
        put(o, "pin", np.zeros(nO)); put(o, "ozmixm", np.zeros((ncell, nO, nM)))
        put(o, "aerosols", np.zeros((ncell, 1, 1))); put(o, "m_hybi", np.zeros((ncell, 1))); put(o, "m_ps", np.zeros(ncell))
        put(o, "cldfrac", np.full((ncell, nV), 0.5)); put(o, "o3clim", np.zeros((ncell, nO)))
        put(o, "o3vmr", np.zeros((ncell, nV)))
        for n in ["rre_cloud", "rre_ice", "rre_snow"]:
            put(o, n, np.full((ncell, nV), 7.0))
        for n in ["glw", "lwcf", "lwdnb", "lwdnbc", "lwdnt", "lwdntc", "lwupb", "lwupbc", "lwupt", "lwuptc", "olrtoa",
                  "coszr", "gsw", "swcf", "swdnb", "swdnbc", "swdnt", "swdntc", "swupb", "swupbc", "swupt", "swuptc",
                  "swddir", "swddni", "swddif"]:
            put(o, n, np.full(ncell, -1.0))
        put(o, "rthratenlw", np.full((ncell, nV), 1.0)); put(o, "rthratensw", np.full((ncell, nV), 1.0))
        for n, a in [("@pres_hyd_p", pres), ("@pres2_hyd_p", pres2), ("@pi_p", pi), ("@t_p", T), ("@t2_p", T2),
                     ("@dz_p", dz), ("@pres_p", pres * rng.uniform(0.999, 1.001, pres.shape)), ("@qv_p", qv),
                     ("@qc_p", qc), ("@qi_p", qi), ("@qs_p", qs)]:
            put(o, n, np.ascontiguousarray(a.T))  # Fortran (nCellsSolve, levels)


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
    bad, total, phase, step = 0, 0, "", 0
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
