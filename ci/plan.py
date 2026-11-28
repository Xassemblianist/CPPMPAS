#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Split the CPU suite into test jobs (.github/workflows/cpu_suite.yml).

usage: plan.py <tests.json (ctest --show-only=json-v1)> <short|long|all> <out dir> [--target SECONDS] [--jobs N]
               [--src-root DIR]

Selection: "short" is every test without the label "long", "long" every test with it, "all" both.
- Tests labelled "gpu" go to one job of their own ("gpu"). It runs them without fixture setups (ctest -FA): this CPU
  build has no CUDA backend, so each one exits 77 and CTest reports it as Skipped (tests/gpu_test.py).
- The other tests are grouped by their fixtures: a test and the setup test of every fixture it requires stay in one
  job, so a model run is made once and its comparisons read its output. The two cheap fixtures "inputs" and
  "work_dir" (make_inputs, make_work_dir) are left out of the grouping; CTest adds their setup test to any job that
  needs it.
- Each group needs the reference releases its tests name (mpas842-dbl-*, ci/data.tsv): in the test's own entry
  (command, arguments, properties) and in the text of every script of the source tree its command runs (a path under
  /w/src, read from --src-root: a script such as tests/check_ref_recipes.py names the releases it reads itself and
  gets only the references directory as its argument). Groups with the same releases
  are packed into jobs of at most --target seconds (estimated wall time at -j PARALLEL: the longest setup plus its
  longest comparison, or the summed time / PARALLEL), from ci/costs.tsv (seconds per test on a GitHub runner;
  unknown tests DEFAULT_COST). A small group (< SMALL seconds) joins the job it adds the fewest bytes of downloads to.

Caching: only the inputs of at most CACHE_MAX_BYTES (compressed) get an Actions cache key ("common" and the small,
often used references: mw, sc); a matrix entry's key<i> is empty for a bigger one. The workflow restores what has a
key, downloads every input the cache did not give from its release (a miss is never a failure), and checks each one
against its SHA-256 (ci/data.sh). The big references are downloaded fresh in every job that needs them: a release
download is about as fast as a cache restore on a GitHub runner, and the cache (10 GB for the repository) then holds
the image and about 1 GB of data instead of every release.

Writes <out>/matrix.json (the jobs: name, tests file, reference tags, timeout), <out>/data.json (the items to fetch),
<out>/<job>.txt (test names, one per line, for ctest --tests-from-file), <out>/expected.tsv (name, labels, job of
every selected test) and <out>/plan.txt (a readable summary).
"""
import argparse
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
LIGHT_FIXTURES = {"inputs", "work_dir"}
DEFAULT_COST = 2.0
PARALLEL = 4
SMALL = 60.0
MAX_REFS = 6  # the workflow's data slots per job
CACHE_MAX_BYTES = 500_000_000  # larger inputs are not cached (module docstring)
TAG_RE = re.compile(r"mpas842-dbl-[a-z0-9]+-[0-9]{3}")
SRC_IN_CONTAINER = "/w/src/"  # ci/in_container.sh mounts the CPPMPAS checkout there


def prop(test, name):
    for p in test.get("properties", []):
        if p["name"] == name:
            return p["value"]
    return None


def load_table(path, columns):
    rows = []
    for line in path.read_text().splitlines():
        if line and not line.startswith("#"):
            parts = line.split("\t")
            if len(parts) != columns:
                sys.exit(f"plan.py: {path.name}: expected {columns} tab-separated columns: {line!r}")
            rows.append(parts)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tests_json")
    ap.add_argument("labels", choices=["short", "long", "all"])
    ap.add_argument("out")
    ap.add_argument("--target", type=float, default=1200.0)
    ap.add_argument("--src-root", type=Path, default=None,
                    help="the checkout mounted at /w/src (its scripts are read for the releases they name)")
    args = ap.parse_args()

    data_rows = load_table(HERE / "data.tsv", 5)
    item_bytes = {}
    for item, _file, _sha, size, _src in data_rows:
        item_bytes[item] = item_bytes.get(item, 0) + int(size)
    costs = {name: float(sec) for name, sec in load_table(HERE / "costs.tsv", 2)}

    tests = json.loads(Path(args.tests_json).read_text())["tests"]
    by_name = {t["name"]: t for t in tests}
    order = {t["name"]: i for i, t in enumerate(tests)}
    if len(by_name) != len(tests):
        sys.exit("plan.py: duplicate test names")
    setup_of = {}
    for t in tests:
        for f in prop(t, "FIXTURES_SETUP") or []:
            setup_of.setdefault(f, []).append(t["name"])

    def labels(name):
        return prop(by_name[name], "LABELS") or []

    def selected(name):
        long_ = "long" in labels(name)
        return args.labels == "all" or (args.labels == "long") == long_

    chosen = [t["name"] for t in tests if selected(t["name"])]
    gpu = [n for n in chosen if "gpu" in labels(n)]
    cpu = [n for n in chosen if "gpu" not in labels(n)]
    if not cpu:
        sys.exit(f"plan.py: labels={args.labels} selects no CPU test")

    # groups: union-find over the fixtures (setups of required fixtures are pulled in, as CTest does)
    parent = {}

    def find(x):
        while parent.setdefault(x, x) != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        parent[find(a)] = find(b)

    members = set(cpu)
    todo = list(cpu)
    while todo:
        n = todo.pop()
        find(n)
        for f in prop(by_name[n], "FIXTURES_REQUIRED") or []:
            if f in LIGHT_FIXTURES:
                continue
            for s in setup_of.get(f, []):
                if "gpu" in labels(s):
                    sys.exit(f"plan.py: CPU test {n} requires the GPU fixture {f}")
                union(n, s)
                if s not in members:
                    members.add(s)
                    todo.append(s)
    groups = {}
    for n in members:
        groups.setdefault(find(n), []).append(n)

    script_tags = {}

    def tags_of_scripts(test):  # the releases named in the scripts of the source tree the command runs
        tags = set()
        for arg in test.get("command", []):
            for word in arg.split():  # also the paths inside a "sh -c" string, quoted there
                word = word.strip("'\";")
                if not word.startswith(SRC_IN_CONTAINER):
                    continue
                if args.src_root is None:
                    sys.exit("plan.py: --src-root is needed (a test runs a script of the source tree)")
                path = args.src_root / word[len(SRC_IN_CONTAINER):]
                if path not in script_tags:
                    script_tags[path] = (set(TAG_RE.findall(path.read_text(errors="replace")))
                                         if path.is_file() else set())
                tags |= script_tags[path]
        return tags

    def refs_of(names):
        tags = set()
        for n in names:
            tags.update(TAG_RE.findall(json.dumps(by_name[n])))
            tags.update(tags_of_scripts(by_name[n]))
        for tag in tags:
            if tag not in item_bytes:
                sys.exit(f"plan.py: {tag} (named by a test) is not in ci/data.tsv: add its archive and SHA-256")
        return frozenset(tags)

    def cost(n):
        return costs.get(n, DEFAULT_COST)

    def is_setup(n):
        return bool(prop(by_name[n], "FIXTURES_SETUP"))

    def chain(names):  # the longest setup and the longest test after it
        setups = [cost(n) for n in names if is_setup(n)]
        rest = [cost(n) for n in names if not is_setup(n)]
        return (max(setups) if setups else 0.0) + (max(rest) if rest else 0.0)

    def wall(names):
        return max(chain(names), sum(cost(n) for n in names) / PARALLEL)

    group_list = [(refs_of(g), sorted(g)) for g in groups.values()]
    group_list.sort(key=lambda rg: (-wall(rg[1]), sorted(rg[0]), rg[1][0]))
    jobs = []  # each: {"refs": set, "tests": [...]}
    small = []
    for refs, names in group_list:
        if wall(names) < SMALL:
            small.append((refs, names))
            continue
        for job in jobs:  # first fit among jobs with the same releases
            if job["refs"] == refs and wall(job["tests"] + names) <= args.target:
                job["tests"] += names
                break
        else:
            jobs.append({"refs": set(refs), "tests": list(names)})
    for refs, names in small:
        best = None
        for job in jobs:
            new = job["refs"] | refs
            if len(new) > MAX_REFS or wall(job["tests"] + names) > args.target:
                continue
            key = (sum(item_bytes[t] for t in new - job["refs"]), wall(job["tests"]))
            if best is None or key < best[0]:
                best = (key, job)
        if best is None:
            jobs.append({"refs": set(refs), "tests": list(names)})
        else:
            best[1]["refs"] |= refs
            best[1]["tests"] += names

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    matrix, names_used, expected = [], set(), []

    def short_tag(tag):
        return tag.replace("mpas842-dbl-", "")

    jobs.sort(key=lambda j: -wall(j["tests"]))
    counter = {}
    for job in jobs:
        if len(job["refs"]) > MAX_REFS:
            sys.exit(f"plan.py: a job needs {len(job['refs'])} releases (at most {MAX_REFS})")
        base = "+".join(short_tag(t) for t in sorted(job["refs"])) or "nodata"
        if len(base) > 40:
            base = f"{len(job['refs'])}refs"
        counter[base] = counter.get(base, 0) + 1
        job["name"] = f"{base}-{counter[base]}"
    if gpu:
        jobs.append({"name": "gpu", "refs": set(), "tests": list(gpu), "gpu": True})
    for job in jobs:
        job["tests"].sort(key=order.__getitem__)
        (out / f"{job['name']}.txt").write_text("".join(n + "\n" for n in job["tests"]))
        est = 30.0 if job.get("gpu") else wall(job["tests"])
        refs = sorted(job["refs"])
        entry = {"name": job["name"], "tests": f"{job['name']}.txt", "gpu": bool(job.get("gpu")),
                 "count": len(job["tests"]), "estimate_s": round(est),
                 "timeout_min": int(min(330, max(60, 3 * est / 60 + 45)))}
        for i in range(MAX_REFS):
            entry[f"ref{i + 1}"] = refs[i] if i < len(refs) else ""
            entry[f"key{i + 1}"] = (cache_key(data_rows, refs[i])
                                    if i < len(refs) and item_bytes[refs[i]] <= CACHE_MAX_BYTES else "")
        matrix.append(entry)
        for n in job["tests"]:
            if n in names_used and n in chosen:
                sys.exit(f"plan.py: {n} is in two jobs")
            names_used.add(n)
            expected.append((n, ";".join(labels(n)), job["name"], "selected" if n in chosen else "setup"))
    missing = set(chosen) - names_used
    if missing:
        sys.exit(f"plan.py: tests in no job: {sorted(missing)[:10]}")

    items = sorted({t for j in jobs for t in j["refs"]})
    data = [{"item": "common", "key": f"cpu-suite-common-{common_key(data_rows)}"}]
    data += [{"item": t, "key": cache_key(data_rows, t) if item_bytes[t] <= CACHE_MAX_BYTES else ""}
             for t in items]
    (out / "matrix.json").write_text(json.dumps({"include": matrix}))
    (out / "data.json").write_text(json.dumps({"include": data}))
    (out / "expected.tsv").write_text("".join("\t".join(r) + "\n" for r in expected))
    lines = [f"labels={args.labels}: {len(chosen)} tests ({len(cpu)} CPU, {len(gpu)} GPU), {len(jobs)} jobs, "
             f"{len(members) - len(cpu)} setup tests added"]
    for m, job in zip(matrix, jobs):
        gb = sum(item_bytes[t] for t in job["refs"]) / 1e9
        refs = [t + ("" if item_bytes[t] <= CACHE_MAX_BYTES else "(fetched)") for t in sorted(job["refs"])]
        lines.append(f"  {m['name']:<28} {m['count']:5d} tests  ~{m['estimate_s']:5d} s  {gb:5.2f} GB  "
                     f"timeout {m['timeout_min']} min  {' '.join(refs)}")
    (out / "plan.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def item_sha(rows, item):
    return next(r[2] for r in rows if r[0] == item)


def cache_key(rows, item):
    return f"cpu-suite-{item}-{item_sha(rows, item)[:16]}"


def common_key(rows):
    import hashlib
    h = hashlib.sha256("".join(f"{r[1]} {r[2]}\n" for r in rows if r[0] == "common").encode())
    return h.hexdigest()[:16]


if __name__ == "__main__":
    main()
