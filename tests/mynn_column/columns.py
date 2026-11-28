# SPDX-License-Identifier: GPL-3.0-or-later
"""Synthetic columns for the MYNN PBL column check, and the comparison of the two programs' outputs.

  columns.py make <in.bin> <nCellsSolve> <nVertLevels> <nsteps> <seed> <restart> <dacycling> <cu> <radt>
                  <scaleaware> <dheat_opt> <len_disp>
  columns.py compare <in.bin> <ref out.bin> <cpp out.bin>   (the input gives nsteps: every step must be in both)

restart: config_do_restart (initflag = 0 from the first step; the MYNN state of in.bin is a previous run's);
dacycling: config_do_DAcycling; cu: 0 'off', 1 'cu_grell_freitas', 2 'cu_kain_fritsch' (physics_get_tend's
convection branches); radt: 0 'off', 1 'rrtmg_lw' / 'rrtmg_sw'; scaleaware, dheat_opt: config_mynn_scaleaware,
config_mynn_dheat_opt (MPAS's defaults are 1); len_disp: config_len_disp (m; dx = len_disp / meshDensity**0.25).

Both files are native-endian streams of records: name (32 bytes, blank padded), int32 kind (0 float64, 1 int32),
int32 ndim, int32 dims[4], then the data in Fortran order. in.bin starts with the records "hdr" (int32: nCellsSolve,
nCells, nEdgesSolve, nEdges, nVertLevels, nsteps, restart, dacycling, cu, radt, scaleaware, dheat_opt) and "dtl"
(float64: dt = 600 s, config_len_disp), then the pool arrays (their Registry names; the harnesses register each in its
pool). The outputs are, per step, the snapshots "A" (after MPAS_to_physics: every pool array and the mpas_atmphys_vars
arrays, names starting with "@"), "C" (after driver_pbl: the plume diagnostics kbl_plume_p, maxwidthbl_p, maxmfbl_p,
zbl_plume_p, then every pool array), "D" (after physics_get_tend with rk_step = 1, dynamics_substep = 1) and "E"
(rk_step = 2), each a marker record then the arrays. The surface-layer outputs MYNN reads (ust, hfx, qfx, ch, qsfc,
rmol, wspd, znt) are inputs, consistent with each column's regime; the MYNN state (qke, el_pbl, tsq, qsq, cov, sh3d,
sm3d, cldfrac_bl, qc_bl, qi_bl) is what a previous step leaves (it matters with restart and at later steps).

nCells = nCellsSolve + 2 (two halo cells: their tendencies are read by tend_toEdges, never computed), nEdges =
3 nCells, nEdgesSolve = nEdges - 3. Cell i (1-based) gets regime (i - 1) % 8:
  0 land, deep daytime convection: hot ground (superadiabatic surface layer: DMP_mf's plumes), a moist mixed layer
    up to 1.5-2.8 km (plume condensation, the mass-flux cloud fraction)
  1 land, stable night, shallow (calm, a surface inversion)
  2 sea, unstable, moist, 8-14 m/s, or 17-31 m/s in cells 9-16, 25-32, ... (DMP_mf's ac_wsp < 1 above 10 m/s)
  3 sea, stable, strong wind
  4 sea ice (xland = 1, cold skin under warmer air)
  5 land, a shallow convective layer under cloud water and ice above the PBL (mym_condensation's rh hacks), negative
    and -0.0 mixing ratios from the transport, |w| > 1 m/s in the PBL (Psig_w < 1, fltv2 = -fltv at Psig_w = 0)
  6 sea, stratocumulus: cloud water at the top of a well-mixed layer, ust < 0.001 in the first step
  7 land, weak heating (fltv near DMP_mf's 0.002 threshold), small wind
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
    """What a complete output of the input's nsteps steps holds: every step's phase markers (A, C, D, E), each followed
    by data records, and at least one value."""
    hdr = [x for n, _, _, x in records(inp) if n == "hdr"]
    if len(hdr) != 1:
        return [f"{inp}: no hdr record"]
    nsteps = int(hdr[0][5])
    if nsteps < 1:
        return [f"{inp}: nsteps = {nsteps}"]
    want = {(ph, s) for ph in ('A', 'C', 'D', 'E') for s in range(1, nsteps + 1)}
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


def make(path, nC, nV, nsteps, seed, restart, dacyc, cu, radt, scaleaware, dheat, len_disp):
    rng = np.random.default_rng(seed)
    nCt = nC + 2
    nE = 3 * nCt
    nEs = nE - 3
    ztop = 25000.0
    eta = np.linspace(0.0, 1.0, nV + 1) ** 1.6  # finer near the ground, ~30 m first layer
    shape_c = (nCt + 1,)
    zg = np.empty((nCt + 1, nV + 1))
    zz = np.empty((nCt + 1, nV))
    names3 = ["rho_zz", "theta_m", "exner", "pressure_base", "pressure_p", "uReconstructZonal",
              "uReconstructMeridional"]
    f = {n: np.empty((nCt + 1, nV)) for n in names3}
    w = np.empty((nCt + 1, nV + 1))
    sc = np.zeros((nCt + 1, nV, 8))
    sfc = {n: np.empty(shape_c) for n in ["skintemp", "xland"]}
    d1 = {n: np.empty(shape_c) for n in ["plrad", "hpbl", "hfx", "qfx", "ust", "wspd", "znt", "ch", "qsfc", "rmol",
                                         "delta", "wstar"]}
    st = {n: np.zeros((nCt + 1, nV)) for n in
          ["el_pbl", "cov", "qke", "qke_adv", "qsq", "tsq", "tke_pbl", "sh3d", "sm3d", "cldfrac_bl", "qc_bl", "qi_bl",
           "edmf_a", "edmf_ent", "edmf_qc", "edmf_qt", "edmf_thl", "edmf_w", "sub_thl", "sub_qv", "det_thl", "det_qv",
           "dqke", "qbuoy", "qdiss", "qshear", "qwt", "kzh", "kzm", "kzq"]}
    surface_pressure = np.empty(shape_c)
    rthlw = np.zeros((nCt + 1, nV))
    rthsw = np.zeros((nCt + 1, nV))
    for i in range(nCt + 1):
        reg = i % 8  # array index i = cell i + 1
        land = reg in (0, 1, 4, 5, 7)
        hs = rng.uniform(0.0, 1200.0) if reg in (0, 1, 5, 7) else 0.0
        zg[i] = hs + eta * (ztop - hs) + 5.0 * np.sin(np.pi * eta) * rng.uniform(-1, 1)
        zm = 0.5 * (zg[i, 1:] + zg[i, :-1])
        zz[i] = (ztop / (ztop - hs)) * (1.0 + 0.003 * rng.uniform(-1, 1, nV))
        z = zm - hs
        tsfc = {0: 302.0, 1: 281.0, 2: 299.0, 3: 284.0, 4: 258.0, 5: 292.0, 6: 289.0, 7: 288.0}[reg] + rng.uniform(-2, 2)
        lapse = np.full(nV, 0.0065)
        zi = {0: rng.uniform(1500, 2800), 1: rng.uniform(60, 200), 2: rng.uniform(600, 1200), 3: rng.uniform(150, 400),
              4: rng.uniform(100, 300), 5: rng.uniform(300, 700), 6: rng.uniform(700, 1100), 7: rng.uniform(300, 600)}[reg]
        if reg in (0, 2, 5, 7):  # well mixed (dry adiabatic) below the PBL top
            lapse[z < zi] = 0.0098
            lapse[(z >= zi) & (z < zi + 200)] = -0.004
        if reg in (1, 3, 4):  # surface inversion
            lapse[z < zi] = -0.012 * rng.uniform(0.5, 1.5)
        if reg == 6:  # mixed layer capped by an inversion
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
        rh0 = {0: 0.7, 1: 0.8, 2: 0.88, 3: 0.9, 4: 0.85, 5: 0.75, 6: 0.97, 7: 0.6}[reg]
        rh = rh0 * rng.uniform(0.9, 1.05, nV)
        if reg == 0:  # moister towards the mixed-layer top: the plumes condense
            rh[z < zi] = np.minimum(0.99, rh0 + 0.25 * z[z < zi] / zi)
        qv = np.clip(rh * qs, 1e-8, 0.03)
        qv[zm > 12000.0] = 3e-6
        qc = np.zeros(nV)
        qi = np.zeros(nV)
        if reg == 6:  # stratocumulus deck below the inversion
            qc[(z > zi - 300) & (z < zi)] = rng.uniform(2e-4, 6e-4)
        if reg in (2, 3, 5):
            qc[(z > 2000) & (z < 3000)] = 1e-4 * rng.uniform(0.2, 1)
            qi[(z > 6000) & (z < 8000)] = 2e-5 * rng.uniform(0.2, 1)
        qr = 1e-5 * rng.uniform(0, 1, nV) * (z < 3000)
        qsn = 1e-5 * rng.uniform(0, 1, nV) * (T < 270)
        qg = np.zeros(nV)
        ni = np.where(qi > 0, rng.uniform(1e3, 1e5, nV), 0.0)
        nr = np.where(qr > 0, rng.uniform(1e2, 1e4, nV), 0.0)
        mix = [qv, qc, qr, qi, qsn, qg, ni, nr]
        for n in range(8):
            sc[i, :, n] = mix[n]
        if reg == 5:  # negative values and -0.0 from the transport
            sc[i, ::5, 1:4] = -1e-9 * rng.uniform(0, 1, (len(range(0, nV, 5)), 3))
            sc[i, 2, 0] = -1e-7
            sc[i, 4, :] = -0.0
            sc[i, 6, 6] = -1.0
        if i in (0, nC - 1):  # ni = -0.0 in the first cell and the last (an odd last cell of 18 or more is the
            sc[i, 3, 6] = -0.0  # reference's scalar maxsd with +0.0 as the destination, PI:397: it stays -0.0)
            sc[i, nV - 2, 6] = -0.0
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
        spd = {0: 3.0, 1: 1.0, 2: 11.0 if i % 16 < 8 else 24.0, 3: 14.0, 4: 6.0, 5: 4.0, 6: 7.0, 7: 1.5}[reg]
        spd = spd * rng.uniform(0.7, 1.3)
        prof = np.log(1.0 + z / 0.1) / np.log(1.0 + 1000.0 / 0.1)
        dirn = rng.uniform(0, 2 * np.pi) + 0.3 * z / 5000.0
        f["uReconstructZonal"][i] = spd * np.minimum(prof, 1.5 + z / 10000.0) * np.cos(dirn)
        f["uReconstructMeridional"][i] = spd * np.minimum(prof, 1.5 + z / 10000.0) * np.sin(dirn)
        w[i] = 0.05 * rng.uniform(-1, 1, nV + 1)
        if reg == 5:
            w[i][: nV // 4] = rng.choice([1.3, 2.5, -1.4]) * rng.uniform(0.8, 1.2)
        w[i, 0] = 0.0
        surface_pressure[i] = ps * rng.uniform(0.999, 1.001)
        dT = {0: 7.0, 1: -5.0, 2: 1.5, 3: -2.5, 4: -4.0, 5: 4.0, 6: 0.8, 7: 0.6}[reg] * rng.uniform(0.6, 1.4)
        sfc["skintemp"][i] = T[0] + 0.0098 * z[0] + dT
        sfc["xland"][i] = 1.0 if land else 2.0
        uns = reg in (0, 2, 5, 6, 7)
        rho0 = rho[0] * (1 + qvp[0])
        ust = {0: 0.45, 1: 0.08, 2: 0.4, 3: 0.55, 4: 0.25, 5: 0.35, 6: 0.0005 if i % 2 else 0.3, 7: 0.15}[reg]
        ust = ust * rng.uniform(0.8, 1.2)
        hfx = {0: 420.0, 1: -25.0, 2: 60.0, 3: -15.0, 4: -40.0, 5: 180.0, 6: 15.0, 7: 3.0}[reg] * rng.uniform(0.7, 1.3)
        qfx = {0: 1.2e-4, 1: 1e-6, 2: 1.5e-4, 3: -2e-6, 4: 2e-6, 5: 6e-5, 6: 5e-5, 7: 1e-5}[reg] * rng.uniform(0.7, 1.3)
        d1["ust"][i] = ust
        d1["hfx"][i] = hfx
        d1["qfx"][i] = qfx
        d1["wspd"][i] = max(spd * 0.5, 0.1)
        d1["znt"][i] = rng.uniform(0.05, 0.8) if land else rng.uniform(1e-4, 5e-4)
        d1["ch"][i] = rng.uniform(0.005, 0.03)
        d1["qsfc"][i] = qv[0] * rng.uniform(1.0, 1.3)
        tv = T[0] * (1 + 0.608 * qvp[0])
        mol = -0.4 * G / tv * (hfx / (rho0 * CP) + 0.608 * T[0] * qfx / rho0) / max(ust, 0.01) ** 3
        d1["rmol"][i] = mol
        d1["hpbl"][i] = zi
        d1["delta"][i] = rng.uniform(-1, 1)
        d1["wstar"][i] = rng.uniform(-1, 1)
        d1["plrad"][i] = rng.uniform(-1, 1)
        # the MYNN state a previous step leaves (restart) or anything (overwritten at initflag = 1)
        inpbl = z < max(zi, 100.0)
        q2 = np.where(inpbl, (2.0 * ust ** 2 + 0.5 * uns) * (1.0 - z / max(zi, 100.0)) + 0.02, 1e-4)
        st["qke"][i] = q2 * rng.uniform(0.8, 1.2, nV)
        st["qke_adv"][i] = st["qke"][i] * rng.uniform(0.9, 1.1, nV)
        st["tke_pbl"][i] = 0.5 * st["qke"][i]
        st["el_pbl"][i] = np.where(inpbl, 0.4 * z / (1 + 0.4 * z / 80.0), 10.0) * rng.uniform(0.9, 1.1, nV)
        st["tsq"][i] = 1e-3 * rng.uniform(0, 1, nV)
        st["qsq"][i] = 1e-8 * rng.uniform(0, 1, nV)
        st["cov"][i] = 1e-6 * rng.uniform(-1, 1, nV)
        st["sh3d"][i] = rng.uniform(0.1, 0.8, nV)
        st["sm3d"][i] = rng.uniform(0.1, 0.8, nV)
        st["cldfrac_bl"][i] = np.where(qc > 0, rng.uniform(0.1, 0.9, nV), 0.0)
        st["qc_bl"][i] = np.where(qc > 0, 1e-5, 0.0)
        st["qi_bl"][i] = np.where(qi > 0, 1e-6, 0.0)
        for n in ["edmf_a", "edmf_ent", "edmf_qc", "edmf_qt", "edmf_thl", "edmf_w", "sub_thl", "sub_qv", "det_thl",
                  "det_qv", "dqke", "qbuoy", "qdiss", "qshear", "qwt", "kzh", "kzm", "kzq"]:
            st[n][i] = rng.uniform(-1, 1, nV)
        if radt:
            rthlw[i] = -2e-5 * rng.uniform(0.5, 1.5, nV)
            rthsw[i] = 1e-5 * rng.uniform(0, 1, nV)
        else:
            rthlw[i] = 1e-5 * rng.uniform(-1, 1, nV)
            rthsw[i] = 1e-5 * rng.uniform(-1, 1, nV)
    coe = rng.integers(1, nCt + 1, size=(nE + 1, 2)).astype(np.int32)
    env = rng.normal(size=(nE + 1, 3))
    env /= np.linalg.norm(env, axis=1, keepdims=True)
    east = rng.normal(size=(nCt + 1, 3))
    east /= np.linalg.norm(east, axis=1, keepdims=True)
    north = rng.normal(size=(nCt + 1, 3))
    north -= (north * east).sum(1, keepdims=True) * east
    north /= np.linalg.norm(north, axis=1, keepdims=True)
    with open(path, "wb") as o:
        put(o, "hdr", np.array([nC, nCt, nEs, nE, nV, nsteps, restart, dacyc, cu, radt, scaleaware, dheat],
                               dtype=np.int32), kind=1)
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
        for n in ["exner", "pressure_base", "pressure_p", "uReconstructZonal", "uReconstructMeridional"]:
            put(o, n, f[n])
        put(o, "rtheta_base", rng.uniform(-1, 1, (nCt + 1, nV)))
        put(o, "rtheta_p", rng.uniform(-1, 1, (nCt + 1, nV)))
        put(o, "rho_edge", rng.uniform(0.1, 1.3, (nE + 1, nV)))
        put(o, "tend_u_phys", rng.uniform(-1, 1, (nE + 1, nV)))
        # sfc_input
        for n in ["skintemp", "xland"]:
            put(o, n, sfc[n])
        # diag_physics
        for n, a in d1.items():
            put(o, n, a)
        put(o, "kpbl", rng.integers(-5, 5, nCt + 1).astype(np.int32), kind=1)
        for n, a in st.items():
            put(o, n, a)
        put(o, "exch_h", rng.uniform(-1, 1, (nCt + 1, nV)))
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
        put(o, "scalars_tend", rng.uniform(-1, 1, (nCt + 1, nV, 8)))
        # the dynamics' arrays of physics_get_tend
        put(o, "tend_ru_physics", rng.uniform(-1, 1, (nE + 1, nV)))
        put(o, "tend_rtheta_physics", rng.uniform(-1, 1, (nCt + 1, nV)))
        put(o, "tend_rho_physics", rng.uniform(-1, 1, (nCt + 1, nV)))


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
        make(sys.argv[2], *map(int, sys.argv[3:14]))
    else:
        sys.exit(compare(sys.argv[2], sys.argv[3], sys.argv[4]))
