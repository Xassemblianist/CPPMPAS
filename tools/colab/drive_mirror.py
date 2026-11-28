#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""The Google Drive mirror of the Xassemblianist/cppmpas-data releases, for GPU runtimes (Colab) that would otherwise
download the release assets from GitHub. Each release asset is kept as it is (split archives stay in their
.partNN assets) at cppwrf-drive:cppmpas-data/<tag>/<asset>; tools/colab/drive_manifest.tsv lists every mirrored
asset with its bytes, SHA-256 (the release's), MD5 (of the bytes uploaded) and Drive path. tools/colab/drive_fetch.sh
reads the mirror back.

  drive_mirror.py plan  [TAG ...]            what a sync would copy (GitHub and Drive metadata only, no transfer)
  drive_mirror.py sync  [--par N] [TAG ...]  copy the missing or stale assets GitHub -> cloud box -> Drive
  drive_mirror.py check [--deep] [TAG ...]   every manifest row against GitHub and Drive; --deep also reads each
                                             Drive copy back on the cloud box and recomputes its SHA-256

TAG: release tags (default: every release except SKIP below; for `inputs` and `mpas-build` only the assets that
ci/data.tsv names). An asset is up to date when the manifest has it with the release's size and SHA-256 and Drive
holds a file of that size and MD5 at its path; sync copies only the others, so running it again after a new
release copies only the new release's assets.

Runs on the workstation, which only reads metadata: the data go from GitHub to the cloud box (~/bin/cppbox,
directory ~/$CPPMPAS_BOX_DIR there, default ~/cppmpas_drive), are checked there against the release's SHA-256 (the asset digest GitHub
reports, and the release's <tag>.sha256 list for split archives), uploaded with rclone, compared with the MD5 Drive
reports (rclone 1.60 gets only MD5 from Drive), and deleted. The box gets only the [cppwrf-drive] section of the
rclone config (scope drive.file), mode 600, removed when the run ends.

The box side is the same file: `drive_mirror.py box-copy` and `box-hash` (internal).
"""
import argparse
import concurrent.futures as cf
import fnmatch
import hashlib
import json
import os
import subprocess
import sys
import threading
import tempfile
import time
import urllib.request

REPO = "Xassemblianist/cppmpas-data"
URL = f"https://github.com/{REPO}/releases/download"
REMOTE = "cppwrf-drive"
ROOT = "cppmpas-data"
HERE = os.path.dirname(os.path.abspath(__file__))
MANIFEST = os.path.join(HERE, "drive_manifest.tsv")
DATA_TSV = os.path.join(HERE, "..", "..", "ci", "data.tsv")
COLS = ("tag", "asset", "bytes", "sha256", "md5", "drive_path")
BOX = os.path.expanduser("~/bin/cppbox")
SCP = os.path.expanduser("~/bin/cppbox-scp")
BOXDIR = os.environ.get("CPPMPAS_BOX_DIR", "cppmpas_drive")  # under the box home
# rclone's shared client ID is rate limited: few calls (256 MiB upload chunks), each retried
RC_FLAGS = ["--low-level-retries", "20", "--drive-chunk-size", "256M"]
# Releases not mirrored (a TAG argument still selects one):
SKIP = (
    "jobs",                                   # CI bundles, not data
    "mpas842-dbl-mw-002",                     # superseded by mw-001 / mw-003, not used by the tests
    "mpas842-dbl-rg-00[1-3]",                 # superseded: made without the libmvec preload (docs/REFERENCES.md)
    "mpas842-dbl-rg-00[5-6]",                 # RG reruns, the tests use rg-004
    "mpas842-dbl-rgcp-00[2-4]",               # determinism / AMD cross-checks of rgcp-001
    "mpas842-dbl-rgcp24-00[2-9]",             # AMD cross-check of rgcp24-001
    "mpas842-dbl-rg60-*", "mpas842-dbl-rg60l-*",    # uniform 60 km / 120 km meshes (about 73 GB): not used by a
    "mpas842-dbl-rg120-*", "mpas842-dbl-rg120l-*",  # test or a GPU run yet; name the tag to mirror one
)
# Releases of which only the assets ci/data.tsv names are mirrored (the rest are WPS / static data no test reads)
PARTIAL = ("inputs", "mpas-build")


def log(*a):
    print("[drive_mirror]", time.strftime("%H:%M:%S"), *a, flush=True)


def run(cmd, timeout, inp=None, check=True):
    r = subprocess.run(cmd, input=inp, capture_output=True, text=True, timeout=timeout)
    if check and r.returncode != 0:
        raise SystemExit(f"[drive_mirror] {' '.join(cmd[:3])} ... failed ({r.returncode}): {r.stderr.strip()[-400:]}")
    return r


def read_manifest():
    m = {}
    if os.path.exists(MANIFEST):
        for line in open(MANIFEST):
            if line.startswith("#") or not line.strip():
                continue
            v = line.rstrip("\n").split("\t")
            if len(v) != len(COLS):
                raise SystemExit(f"[drive_mirror] bad manifest line: {line!r}")
            m[(v[0], v[1])] = dict(zip(COLS, v))
    return m


def write_manifest(m):
    rows = sorted(m.values(), key=lambda r: (r["tag"], r["asset"]))
    with open(MANIFEST + ".new", "w") as f:
        f.write("# The Google Drive mirror of the Xassemblianist/cppmpas-data release assets (tools/colab/drive_mirror.py\n"
                "# writes it, tools/colab/drive_fetch.sh reads it). sha256: the release asset's; md5: Drive's, of the\n"
                "# same bytes; drive_path: in the rclone remote cppwrf-drive (scope drive.file).\n")
        f.write("# " + "\t".join(COLS) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in COLS) + "\n")
    os.replace(MANIFEST + ".new", MANIFEST)


def releases():
    r = run(["gh", "api", "--paginate", f"repos/{REPO}/releases?per_page=100",
             "--jq", ".[] | {tag: .tag_name, assets: [.assets[] | {name, size, digest}]}"], 300)
    return [json.loads(l) for l in r.stdout.splitlines() if l.strip()]


def partial_names():
    names = {}
    for line in open(DATA_TSV):
        if line.startswith("#") or not line.strip():
            continue
        v = line.rstrip("\n").split("\t")
        if v[4].startswith("release:"):
            names.setdefault(v[4][len("release:"):], set()).add(v[1])
    return names


def wanted(tags):
    """[(tag, asset, bytes, sha256)] of the selected releases, each SHA-256 the asset digest GitHub reports, and for
    split archives also the line of the release's <tag>.sha256 list (both must agree)."""
    rel = {r["tag"]: r for r in releases()}
    if tags:
        miss = [t for t in tags if t not in rel]
        if miss:
            raise SystemExit(f"[drive_mirror] no release {' '.join(miss)}")
        sel = tags
    else:
        sel = [t for t in rel if not any(fnmatch.fnmatchcase(t, p) for p in SKIP)]
    part = partial_names()
    out = []
    for t in sorted(sel):
        assets = rel[t]["assets"]
        if t in PARTIAL:
            assets = [a for a in assets if a["name"] in part.get(t, ())]
            if len(assets) != len(part.get(t, ())):
                raise SystemExit(f"[drive_mirror] {t}: ci/data.tsv names assets the release lacks")
        lists = {}
        for a in assets:
            if a["name"].endswith(".sha256"):  # the release's own lists (a few hundred bytes)
                with urllib.request.urlopen(f"{URL}/{t}/{a['name']}", timeout=60) as f:
                    for line in f.read().decode().splitlines():
                        if line.strip():
                            d, n = line.split()
                            lists[n.lstrip("*")] = d.lower()
        for a in assets:
            dg = a.get("digest") or ""
            sha = dg[len("sha256:"):] if dg.startswith("sha256:") else lists.get(a["name"], "")
            if not sha:
                raise SystemExit(f"[drive_mirror] {t}/{a['name']}: no SHA-256 (no GitHub digest, not in a .sha256 list)")
            if a["name"] in lists and lists[a["name"]] != sha:
                raise SystemExit(f"[drive_mirror] {t}/{a['name']}: GitHub digest {sha} differs from the release's list")
            out.append((t, a["name"], int(a["size"]), sha))
    return out


def drive_listing():
    run(["rclone", "mkdir", f"{REMOTE}:{ROOT}"], 120)  # rclone 1.60 can hang listing a missing folder
    r = run(["rclone", "lsjson", "-R", "--files-only", "--hash", f"{REMOTE}:{ROOT}"], 600)
    d = {}
    for e in json.loads(r.stdout):
        p = f"{ROOT}/{e['Path']}"
        # Drive allows two files of one name: such a path matches nothing
        d[p] = (-1, "duplicate") if p in d else (int(e["Size"]), e.get("Hashes", {}).get("md5", ""))
    return d


def status(want, man, drive):
    todo, ok = [], []
    for t, a, b, sha in want:
        p = f"{ROOT}/{t}/{a}"
        r = man.get((t, a))
        good = (r is not None and int(r["bytes"]) == b and r["sha256"] == sha and r["drive_path"] == p and
                drive.get(p) == (b, r["md5"]))
        (ok if good else todo).append((t, a, b, sha))
    return todo, ok


def box(cmd, timeout):
    return subprocess.run([BOX, cmd], capture_output=True, text=True, timeout=timeout)


def cmd_plan(args):
    want = wanted(args.tags)
    todo, ok = status(want, read_manifest(), drive_listing())
    for t, a, b, _ in todo:
        print(f"COPY {t}/{a} {b}")
    log(f"{len(want)} assets: {len(ok)} up to date, {len(todo)} to copy ({sum(x[2] for x in todo) / 1e9:.2f} GB)")
    return 0


def with_box_conf(fn):
    """Run fn() with the [cppwrf-drive] rclone section and this script in ~/BOXDIR on the box; remove the section
    afterwards (also when fn fails or is interrupted)."""
    r = box(f"umask 077; mkdir -p ~/{BOXDIR} && chmod 700 ~/{BOXDIR} && "
            f"! {{ [ -f ~/{BOXDIR}/job.pid ] && kill -0 $(cat ~/{BOXDIR}/job.pid) 2> /dev/null; }} && echo READY", 300)
    if "READY" not in r.stdout:
        raise SystemExit("[drive_mirror] the box is not reachable or a job is still running in ~/" + BOXDIR)
    tmp = os.path.join(tempfile.mkdtemp(), f"rclone.{os.getpid()}.conf")  # mode 700 directory
    try:
        old = os.umask(0o077)
        try:
            with open(tmp, "w") as f:
                subprocess.run(["rclone", "config", "show", REMOTE], stdout=f, check=True, timeout=60)
        finally:
            os.umask(old)
        subprocess.run([SCP, tmp, os.path.abspath(__file__), f"cppbox:{BOXDIR}/"], check=True, timeout=300,
                       capture_output=True)
        os.remove(tmp)
        box(f"cd ~/{BOXDIR} && mv {os.path.basename(tmp)} rclone.conf && chmod 600 rclone.conf", 120)
        return fn()
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)
        os.rmdir(os.path.dirname(tmp))
        r = box(f"rm -f ~/{BOXDIR}/rclone.conf ~/{BOXDIR}/rclone.*.conf; test ! -e ~/{BOXDIR}/rclone.conf && echo GONE", 300)
        if "GONE" in r.stdout:
            log("rclone section removed from the box")
        else:
            log(f"WARNING: could not confirm the removal of ~/{BOXDIR}/rclone.conf on the box")


def box_job(sub, todo, extra, poll):
    """Start `drive_mirror.py SUB` detached on the box on the rows `todo`, follow its log, return its result rows."""
    lines = "".join("\t".join(map(str, x)) + "\n" for x in todo)
    r = subprocess.run([BOX, f"cat > ~/{BOXDIR}/todo.tsv"], input=lines, text=True, capture_output=True, timeout=300)
    if r.returncode:
        raise SystemExit("[drive_mirror] could not write the box's todo list")
    # only the job goes to the background (with its own output), so the ssh session returns at once
    box(f"cd ~/{BOXDIR} || exit 1; rm -f results.tsv job.log; setsid nohup python3 -u drive_mirror.py {sub} "
        f"--conf rclone.conf --todo todo.tsv --out results.tsv {extra} < /dev/null > job.log 2>&1 & echo $! > job.pid", 120)
    seen = 0
    while True:
        time.sleep(poll)
        r = box(f"cd ~/{BOXDIR} && cat job.log; kill -0 $(cat job.pid) 2> /dev/null && echo JOB-RUNNING || echo JOB-ENDED",
                180)
        out = [l for l in r.stdout.splitlines() if l.startswith("[box]")]
        for l in out[seen:]:
            print(l, flush=True)
        seen = max(seen, len(out))
        if "JOB-ENDED" in r.stdout:
            break
    r = box(f"cat ~/{BOXDIR}/results.tsv", 180)
    return [l.split("\t") for l in r.stdout.splitlines() if l.strip()]


def cmd_sync(args):
    want = wanted(args.tags)
    man = read_manifest()
    todo, ok = status(want, man, drive_listing())
    log(f"{len(want)} assets: {len(ok)} up to date, {len(todo)} to copy ({sum(x[2] for x in todo) / 1e9:.2f} GB)")
    if todo:
        todo.sort(key=lambda x: -x[2])
        res = with_box_conf(lambda: box_job("box-copy", todo, f"--par {args.par}", args.poll))
        for v in res:
            if v[0] == "OK":
                t, a, b, sha, md5 = v[1:6]
                man[(t, a)] = dict(tag=t, asset=a, bytes=b, sha256=sha, md5=md5, drive_path=f"{ROOT}/{t}/{a}")
        write_manifest(man)
    todo, ok = status(want, read_manifest(), drive_listing())
    for t, a, _, _ in todo:
        log(f"NOT ON DRIVE: {t}/{a}")
    log(f"{len(ok)} of {len(want)} assets on Drive with the release's size and SHA-256 and the recorded MD5")
    return 1 if todo else 0


def cmd_check(args):
    """Every manifest row (of TAGs): the release still has the asset with that size and SHA-256 and Drive holds a file
    of that size and MD5 at drive_path; with --deep, the box reads each Drive copy back and recomputes its SHA-256."""
    man = read_manifest()
    rows = [r for r in man.values() if not args.tags or r["tag"] in args.tags]
    rel = {(t, a): (b, sha) for t, a, b, sha in wanted(sorted({r["tag"] for r in rows}))}
    drive = drive_listing()
    bad = 0
    res = {}
    if args.deep:
        todo = [(r["tag"], r["asset"], r["bytes"], r["sha256"]) for r in rows]
        todo.sort(key=lambda x: -int(x[2]))
        for v in with_box_conf(lambda: box_job("box-hash", todo, f"--par {args.par}", args.poll)):
            res[(v[1], v[2])] = v
    print("tag\tasset\tbytes\tsha256\tdrive_md5\tverified")
    for r in sorted(rows, key=lambda r: (r["tag"], r["asset"])):
        k = (r["tag"], r["asset"])
        why = []
        if rel.get(k) != (int(r["bytes"]), r["sha256"]):
            why.append("release differs")
        if drive.get(r["drive_path"]) != (int(r["bytes"]), r["md5"]):
            why.append(f"Drive has {drive.get(r['drive_path'])}")
        if args.deep:
            v = res.get(k)
            if not v or v[0] != "OK":
                why.append("re-read SHA-256 " + (v[0] if v else "missing"))
        how = "size+md5" + ("+sha256 re-read" if args.deep else "")
        print("\t".join([r["tag"], r["asset"], r["bytes"], r["sha256"], drive.get(r["drive_path"], ("", "-"))[1],
                         "yes (" + how + ")" if not why else "NO: " + "; ".join(why)]))
        bad += bool(why)
    log(f"{len(rows) - bad} of {len(rows)} manifest assets verified" + (" (Drive copies re-read)" if args.deep else ""))
    return 1 if bad else 0


# ---- box side ------------------------------------------------------------------------------------------------------
def hashes(stream):
    s, m, n = hashlib.sha256(), hashlib.md5(), 0
    while True:
        c = stream.read(1 << 22)
        if not c:
            break
        s.update(c)
        m.update(c)
        n += len(c)
    return n, s.hexdigest(), m.hexdigest()


def box_rows(path):
    return [l.rstrip("\n").split("\t") for l in open(path) if l.strip()]


def drive_stat(rc, t, a):
    """[size, md5] of cppmpas-data/<t>/<a>, or None; from the folder's listing (lsjson on a file path is not reliable
    in rclone 1.60)."""
    r = subprocess.run(rc + ["lsjson", "--hash", "--files-only", f"{REMOTE}:{ROOT}/{t}"], capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError(f"rclone lsjson {ROOT}/{t}: {r.stderr.strip()[-200:]}")
    e = [x for x in json.loads(r.stdout) if x["Name"] == a]
    if len(e) > 1:
        raise RuntimeError(f"{len(e)} files named {ROOT}/{t}/{a} on Drive")
    return [int(e[0]["Size"]), e[0].get("Hashes", {}).get("md5", "")] if e else None


def box_run(args, one):
    lock = threading.Lock()
    rows = box_rows(args.todo)
    rc = ["rclone", "--config", os.path.abspath(args.conf)]
    # Drive allows two folders of one name: parallel uploads into a missing folder each make one. Make them first.
    subprocess.run(rc + ["mkdir", f"{REMOTE}:{ROOT}"], check=True)
    for t in sorted({r[0] for r in rows}):
        subprocess.run(rc + ["mkdir", f"{REMOTE}:{ROOT}/{t}"], check=True)
    with open(args.out, "w") as out, cf.ThreadPoolExecutor(args.par) as ex:
        def job(row):
            t, a, b, sha = row
            for attempt in range(3):
                try:
                    res = one(t, a, int(b), sha)
                    break
                except Exception as e:  # a failed transfer is retried twice, then reported
                    res = ["FAIL", t, a, b, sha, str(e)[-300:]]
                    print(f"[box] {time.strftime('%H:%M:%S')} RETRY {t}/{a}: {e}", flush=True)
            with lock:
                out.write("\t".join(res) + "\n")
                out.flush()
                print(f"[box] {time.strftime('%H:%M:%S')} {res[0]} {t}/{a} {b}", flush=True)
            return res[0] == "OK"
        good = sum(ex.map(job, rows))
    print(f"[box] {time.strftime('%H:%M:%S')} done: {good} of {len(rows)} OK", flush=True)
    return 0 if good == len(rows) else 1


def cmd_box_copy(args):
    work = os.path.abspath("work")
    rc = ["rclone", "--config", os.path.abspath(args.conf)]

    def one(t, a, b, sha):
        d = os.path.join(work, t)
        os.makedirs(d, exist_ok=True)
        f = os.path.join(d, a)
        try:
            t0 = time.time()
            subprocess.run(["curl", "-fsSL", "--retry", "5", "--connect-timeout", "30", "--max-time", "3600", "-o", f,
                            f"{URL}/{t}/{a}"], check=True)
            with open(f, "rb") as s:
                n, s256, md5 = hashes(s)
            if n != b or s256 != sha:
                raise RuntimeError(f"download is {n} B sha256 {s256}, the release says {b} B {sha}")
            t1 = time.time()
            subprocess.run(rc + ["copyto", *RC_FLAGS, "--retries", "3", f, f"{REMOTE}:{ROOT}/{t}/{a}"], check=True,
                           capture_output=True)
            st = drive_stat(rc, t, a)
            if st != [b, md5]:
                raise RuntimeError(f"Drive reports {st}, the upload was {b} B md5 {md5}")
            print(f"[box] {t}/{a}: GitHub {t1 - t0:.0f} s, Drive {time.time() - t1:.0f} s", flush=True)
            return ["OK", t, a, str(b), sha, md5]
        finally:
            if os.path.exists(f):
                os.remove(f)
    return box_run(args, one)


def cmd_box_hash(args):
    rc = ["rclone", "--config", os.path.abspath(args.conf)]

    def one(t, a, b, sha):
        p = subprocess.Popen(rc + ["cat", *RC_FLAGS, f"{REMOTE}:{ROOT}/{t}/{a}"], stdout=subprocess.PIPE)
        n, s256, md5 = hashes(p.stdout)
        if p.wait() != 0:
            raise RuntimeError("rclone cat failed")
        return ["OK" if (n, s256) == (b, sha) else "MISMATCH", t, a, str(n), s256, md5]
    return box_run(args, one)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    for name in ("plan", "sync", "check"):
        p = sp.add_parser(name)
        p.add_argument("tags", nargs="*")
        p.add_argument("--par", type=int, default=2)
        p.add_argument("--poll", type=int, default=60)
        p.add_argument("--deep", action="store_true")
    for name in ("box-copy", "box-hash"):
        p = sp.add_parser(name)
        p.add_argument("--conf", required=True)
        p.add_argument("--todo", required=True)
        p.add_argument("--out", required=True)
        p.add_argument("--par", type=int, default=2)
    a = ap.parse_args()
    return {"plan": cmd_plan, "sync": cmd_sync, "check": cmd_check,
            "box-copy": cmd_box_copy, "box-hash": cmd_box_hash}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
