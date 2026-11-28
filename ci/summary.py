#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""The CPU suite's verdict (.github/workflows/cpu_suite.yml): every selected test of the plan (ci/plan.py) against the
CTest JUnit files of all test jobs.

usage: summary.py <plan dir> <results dir (one sub-directory per job, each with junit.xml)> <out dir>

- A test without the label "gpu" passes only if it ran and passed in every job that ran it (a fixture setup such as
  make_inputs runs in several jobs). Failed, Not Run (a fixture failed, an executable is missing), Skipped, or no
  result at all: the suite fails. A check that did not run is a failure.
- A test with the label "gpu" must be reported Skipped (exit status 77: this CPU build has no CUDA backend). It is
  counted as skipped, never as passed; a GPU test that passes or fails here is an error.
Writes <out>/summary.md (also the job summary), summary.json, failed.txt, and costs_measured.tsv (the longest
measured time per test, the format of ci/costs.tsv). Exit status 1 when the suite fails.
"""
import json
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

SKIP_MESSAGE = "SKIP_RETURN_CODE=77"


def main():
    if len(sys.argv) != 4:
        print(__doc__)
        return 2
    plan, results, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    out.mkdir(parents=True, exist_ok=True)
    expected = {}
    for line in (plan / "expected.tsv").read_text().splitlines():
        name, labels, job, kind = line.split("\t")
        expected[name] = {"labels": labels.split(";") if labels else [], "job": job, "kind": kind}
    jobs = json.loads((plan / "matrix.json").read_text())["include"]

    runs = {}  # name -> [(job, status, seconds, message)]
    job_state = {}
    for job in jobs:
        junit = results / job["name"] / "junit.xml"
        if not junit.is_file():
            job_state[job["name"]] = "no junit.xml (the job failed before or during ctest)"
            continue
        root = ET.parse(junit).getroot()
        cases = root.iter("testcase")
        n = 0
        for case in cases:
            n += 1
            status = case.get("status")
            skipped = case.find("skipped")
            message = skipped.get("message", "") if skipped is not None else ""
            if status == "run":
                st = "passed"
            elif status == "fail":
                st = "failed"
            elif status == "notrun" and message == SKIP_MESSAGE:
                st = "skipped"
            else:
                st = "not run"
            runs.setdefault(case.get("name"), []).append((job["name"], st, float(case.get("time") or 0), message))
        job_state[job["name"]] = f"{n} results"

    verdict = {}
    problems = []
    for name, info in expected.items():
        gpu = "gpu" in info["labels"]
        got = runs.get(name, [])
        states = {st for _, st, _, _ in got}
        if not got:
            v = "missing"
        elif gpu:
            v = "skipped" if states == {"skipped"} else ("ERROR: " + ",".join(sorted(states)))
        elif states == {"passed"}:
            v = "passed"
        elif "failed" in states:
            v = "failed"
        elif "not run" in states:
            v = "not run"
        else:
            v = "skipped (a CPU test that did not run)"
        verdict[name] = v
        ok = v == ("skipped" if gpu else "passed")
        if not ok:
            where = "; ".join(f"{j}: {st} {m}".strip() for j, st, _, m in got) or "no result in any job"
            problems.append(f"{name} [{info['job']}]: {v} ({where})")
    # tests CTest added to a job by itself (fixture setups) must pass as well
    for name, got in runs.items():
        if name not in expected and any(st != "passed" for _, st, _, _ in got):
            problems.append(f"{name} (fixture setup added by CTest): " + "; ".join(f"{j}: {st}" for j, st, _, _ in got))

    def category(info):
        return ("long" if "long" in info["labels"] else "short") + (" GPU" if "gpu" in info["labels"] else " CPU")

    table = {}
    for name, info in expected.items():
        if info["kind"] != "selected":
            continue
        row = table.setdefault(category(info), {"tests": 0, "passed": 0, "skipped": 0, "failed": 0, "not run": 0,
                                                "missing": 0, "error": 0})
        row["tests"] += 1
        v = verdict[name]
        key = ("passed" if v == "passed" else "skipped" if v == "skipped" else "failed" if v == "failed"
               else "not run" if v.startswith("not run") or v.startswith("skipped (") else "missing"
               if v == "missing" else "error")
        row[key] += 1

    ok = not problems
    lines = ["# CPPMPAS CPU suite", "", f"**{'PASSED' if ok else 'FAILED'}**", ""]
    commit = plan / "commit.txt"
    if commit.is_file():
        lines += [f"Commit: `{commit.read_text().strip()}`", ""]
    lines += ["| label | tests | passed | skipped | failed | not run | missing | error |",
              "|---|---:|---:|---:|---:|---:|---:|---:|"]
    for cat in sorted(table):
        r = table[cat]
        lines.append(f"| {cat} | {r['tests']} | {r['passed']} | {r['skipped']} | {r['failed']} | {r['not run']} | "
                     f"{r['missing']} | {r['error']} |")
    lines += ["", "GPU tests: this is a CPU build (no CUDA backend); each GPU test exits 77 and is reported as "
              "skipped, not passed.", "", "| job | results |", "|---|---|"]
    lines += [f"| {j} | {s} |" for j, s in job_state.items()]
    if problems:
        lines += ["", f"## {len(problems)} problem(s)", ""] + [f"- {p}" for p in problems[:300]]
    (out / "summary.md").write_text("\n".join(lines) + "\n")
    (out / "failed.txt").write_text("".join(p + "\n" for p in problems))
    (out / "summary.json").write_text(json.dumps({"passed": ok, "counts": table, "jobs": job_state,
                                                  "problems": problems}, indent=1))
    times = {}
    for name, got in runs.items():
        times[name] = max(t for _, _, t, _ in got)
    (out / "costs_measured.tsv").write_text(
        "".join(f"{n}\t{t:.0f}\n" for n, t in sorted(times.items(), key=lambda kv: -kv[1]) if t >= 5))
    print("\n".join(lines))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
