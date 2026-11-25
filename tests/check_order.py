#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Check CPPMPAS's in-memory mesh (cppmpas dump) against an independent re-implementation of MPAS's one-rank
numbering, written differently on purpose (numpy, no loops shared with the C++):

* local edge / vertex order = order of first appearance in edgesOnCell / verticesOnCell walked cell by cell,
  slots 1..nEdgesOnCell (mpas_block_decomp.F:633-685); cells = file order (mpas_block_decomp.F:303-315);
* connectivity: every slot of elements 1..n holds the local index of the global ID, n+1 when absent
  (mpas_stream_manager.F:5619-5636); the garbage column keeps the Registry default 0;
* other arrays: local element i holds file record indexTo*ID(i); garbage element = Registry default 0.

usage: check_order.py <grid or init file> <dump.nc> [MW|SC|JW]
"""
import sys

import netCDF4 as nc
import numpy as np

# sequences quoted from ncdump of the grid files in NEONM02_verdicts.md F2 / DYCORE_MAP.md §4.12
QUOTED = {
    "MW": {"edges": [1, 2, 3, 4, 14, 24, 5, 6, 7, 17, 15], "vertices": [1, 2, 3193, 4, 3, 10]},
    "SC": {"edges": [1, 2, 3, 469, 5, 83775, 4], "vertices": [1, 2, 623, 314, 313, 4]},
}
CONN = {"cellsOnCell": "c", "edgesOnCell": "e", "verticesOnCell": "v", "cellsOnEdge": "c", "verticesOnEdge": "v",
        "edgesOnEdge": "e", "cellsOnVertex": "c", "edgesOnVertex": "e"}
OUTER = {"nCells": "c", "nEdges": "e", "nVertices": "v"}


def first_seen(on_cell, n_on_cell):
    used = np.arange(on_cell.shape[1])[None, :] < n_on_cell[:, None]
    flat = on_cell[used]                      # C order = cell by cell, slot by slot
    ids, first = np.unique(flat, return_index=True)
    return ids[np.argsort(first)]


def main():
    src_path, dump_path = sys.argv[1], sys.argv[2]
    case = sys.argv[3] if len(sys.argv) > 3 else None
    failures = []
    with nc.Dataset(src_path) as s, nc.Dataset(dump_path) as d:
        s.set_auto_mask(False)
        d.set_auto_mask(False)
        n = {"c": len(s.dimensions["nCells"]), "e": len(s.dimensions["nEdges"]), "v": len(s.dimensions["nVertices"])}
        nEdgesOnCell = s["nEdgesOnCell"][:]
        order = {"c": np.arange(1, n["c"] + 1),
                 "e": first_seen(s["edgesOnCell"][:], nEdgesOnCell),
                 "v": first_seen(s["verticesOnCell"][:], nEdgesOnCell)}
        for k in "ev":
            if order[k].size != n[k]:
                failures.append(f"independent order of {k} has {order[k].size} elements, file {n[k]}")
        inv = {}
        for k, o in order.items():
            i = np.zeros(n[k] + 1, np.int64)
            i[o] = np.arange(1, n[k] + 1)
            inv[k] = i
        for name, k in (("indexToCellID", "c"), ("indexToEdgeID", "e"), ("indexToVertexID", "v")):
            got = d[name][:]
            if not np.array_equal(got[:-1], order[k]) or got[-1] != 0:
                bad = np.flatnonzero(got[:-1] != order[k])
                failures.append(f"{name}: differs from the first-seen order at local {bad[:5] + 1}, garbage {got[-1]}")
        if case in QUOTED:
            for label, k in (("edges", "e"), ("vertices", "v")):
                q = QUOTED[case][label]
                got = list(d["indexToEdgeID" if k == "e" else "indexToVertexID"][: len(q)])
                if got != q:
                    failures.append(f"{case} {label}: first local IDs {got}, quoted {q}")
        checked = 0
        for name in d.variables:
            if name.startswith("indexTo"):
                continue
            if name not in s.variables:
                failures.append(f"{name}: in the dump but not in the source")
                continue
            mem = d[name][:]
            sv = s[name]
            src = sv[0] if sv.dimensions and sv.dimensions[0] == "Time" else sv[:]
            sdims = [x for x in sv.dimensions if x != "Time"]
            outer = OUTER.get(sdims[0]) if sdims else None
            if outer is None:
                if not np.array_equal(mem.reshape(src.shape), src):
                    failures.append(f"{name}: values differ")
                checked += 1
                continue
            body, garbage = mem[:-1], mem[-1]
            expect = src[order[outer] - 1]
            if name in CONN:
                k = CONN[name]
                valid = (expect >= 1) & (expect <= n[k])
                expect = np.where(valid, inv[k][np.clip(expect, 0, n[k])], n[k] + 1)
            if not np.array_equal(body, expect):
                bad = np.argwhere(body != expect)
                failures.append(f"{name}: {len(bad)} values differ, first local index {bad[0] + 1}")
            if np.any(garbage != 0):
                failures.append(f"{name}: garbage element is not the Registry default 0")
            checked += 1
        if checked < 20:
            failures.append(f"only {checked} arrays checked")
    for f in failures:
        print("FAIL", f)
    if failures:
        return 1
    print(f"OK {checked} arrays + 3 index maps, order and connectivity as MPAS on one rank")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
