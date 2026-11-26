#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""tests/mutation_check.py's reading of CTest output: a test that failed or crashed kills the mutant; a test that
passed, was not run because its fixture setup failed (CTest's "***Not Run"), was skipped or timed out does not; nor
does one that did not pass before the mutant was applied. Its main() on a simulated cmake / ctest: the selection is
parsed strictly, CTest runs one test at a time, test names are matched literally, and a test that fails on the
restored source makes the run fail.

usage: test_mutation_check.py
"""
import contextlib
import io
import os
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import mutation_check  # noqa: E402

# a CTEST_PARALLEL_LEVEL in the caller's environment must not make the harness run tests in parallel
os.environ["CTEST_PARALLEL_LEVEL"] = "8"

# CTest 3.31 output of `ctest -R ^dep$` where dep requires a fixture whose setup test failed
NOT_RUN = """Test project /w/b
    Start 1: setup
1/2 Test #1: setup ............................***Failed    0.02 sec
    Start 2: dep
Failed test dependencies: setup
2/2 Test #2: dep ..............................***Not Run   0.00 sec

0% tests passed, 2 tests failed out of 2

The following tests FAILED:
\t  1 - setup (Failed)
\t  2 - dep (Not Run)
Errors while running CTest
"""
PASSED = """    Start 205: initdump_run_JW
1/2 Test #205: initdump_run_JW .......................   Passed    1.55 sec
    Start 207: initdump_JW_ke
2/2 Test #207: initdump_JW_ke ........................   Passed    0.30 sec

100% tests passed, 0 tests failed out of 2
"""
FAILED = """1/2 Test #205: initdump_run_JW .......................   Passed    1.55 sec
2/2 Test #207: initdump_JW_ke ........................***Failed    0.30 sec

The following tests FAILED:
\t207 - initdump_JW_ke (Failed)
"""
SEGV = """1/1 Test  #12: unit_solve_diag .......................***Exception: SegFault  0.10 sec
"""
ABORT = """1/1 Test  #12: unit_solve_diag .......................Subprocess aborted***Exception:   0.10 sec
"""
DISABLED = """1/1 Test  #12: unit_solve_diag .......................***Not Run (Disabled)   0.00 sec
"""


class Status(unittest.TestCase):
    def k(self, out, t):
        return mutation_check.killed(mutation_check.test_status(out, t))

    def test_not_run_is_not_a_kill(self):
        self.assertEqual(mutation_check.test_status(NOT_RUN, "dep"), "***Not Run")
        self.assertFalse(self.k(NOT_RUN, "dep"))
        self.assertFalse(self.k(DISABLED, "unit_solve_diag"))

    def test_failed_setup_itself_is_a_kill(self):
        self.assertTrue(self.k(NOT_RUN, "setup"))

    def test_passed(self):
        self.assertFalse(self.k(PASSED, "initdump_JW_ke"))
        self.assertFalse(self.k(PASSED, "initdump_run_JW"))

    def test_failed_and_crashed(self):
        self.assertTrue(self.k(FAILED, "initdump_JW_ke"))
        self.assertTrue(self.k(SEGV, "unit_solve_diag"))
        self.assertTrue(self.k(ABORT, "unit_solve_diag"))

    def test_baseline(self):
        # a test that fails, crashes or is not run before any mutant is applied cannot show a kill
        self.assertEqual(mutation_check.not_passing(FAILED, ["initdump_run_JW", "initdump_JW_ke"]), {"initdump_JW_ke"})
        self.assertEqual(mutation_check.not_passing(NOT_RUN, ["setup", "dep"]), {"setup", "dep"})
        self.assertEqual(mutation_check.not_passing(PASSED, ["initdump_run_JW", "initdump_JW_ke", "absent"]), {"absent"})

    def test_absent_or_prefix(self):
        self.assertFalse(self.k(FAILED, "initdump_JW"))  # a prefix of a test name is not that test
        self.assertFalse(self.k("", "initdump_JW_ke"))

    def test_skipped_and_timeout_are_not_kills(self):
        self.assertEqual(mutation_check.test_status(SKIPPED, "unit_namelist_locale"), "***Skipped")
        self.assertFalse(self.k(SKIPPED, "unit_namelist_locale"))
        self.assertEqual(mutation_check.test_status(TIMEOUT, "unit_libm"), "***Timeout")
        self.assertFalse(self.k(TIMEOUT, "unit_libm"))

    def test_is_kill(self):
        self.assertTrue(mutation_check.is_kill(8, "***Failed", "t", set()))
        self.assertFalse(mutation_check.is_kill(8, "***Failed", "t", {"t"}))  # failed before the mutant
        self.assertFalse(mutation_check.is_kill(0, "***Failed", "t", set()))  # CTest itself saw no failure
        self.assertFalse(mutation_check.is_kill(8, "Passed", "t", set()))


SKIPPED = """1/1 Test  #87: unit_namelist_locale .............***Skipped   0.01 sec
"""
TIMEOUT = """1/1 Test   #5: unit_libm ........................***Timeout 1500.02 sec
"""


class Selection(unittest.TestCase):
    def test_accepted(self):
        self.assertEqual(mutation_check.parse_selection("0", 5), [0, 1, 2, 3, 4])
        self.assertEqual(mutation_check.parse_selection("4", 5), [4])
        self.assertEqual(mutation_check.parse_selection("3,1", 5), [3, 1])
        self.assertEqual(mutation_check.parse_selection("0,4", 5), [0, 4])

    def test_refused(self):
        for sel in ("122,123,", ",1", "1,,2", "1,", "", "-1", "+1", " 1", "1 ,2", "1;2", "a", "1,1", "5", "1,5"):
            with self.subTest(sel=sel):
                with self.assertRaises(SystemExit) as e:
                    mutation_check.parse_selection(sel, 5)
                self.assertIsInstance(e.exception.code, str)  # a message, not a traceback or an exit status


class Main(unittest.TestCase):
    """main() on a temporary source tree; cmake always builds, ctest is simulated: fails(test, source text, call
    number, parallel) says whether a test fails"""

    def run_main(self, mutations, tests, fails, sel="0"):
        with tempfile.TemporaryDirectory() as d:
            src = Path(d)
            (src / "f.cpp").write_text("int x = 1;\n")
            calls = []

            def run(cmd, **kw):
                if cmd[0] == "cmake":
                    return subprocess.CompletedProcess(cmd, 0, "", "")
                calls.append(cmd)
                regex = cmd[cmd.index("-R") + 1]
                if "-j" in cmd:
                    parallel = int(cmd[cmd.index("-j") + 1]) > 1
                else:
                    parallel = int(os.environ.get("CTEST_PARALLEL_LEVEL", "1")) > 1
                names = [t for t in tests if re.search(regex, t)]
                body = (src / "f.cpp").read_text()
                out, rc = [], 0
                for k, t in enumerate(names, 1):
                    bad = fails(t, body, len(calls), parallel and len(names) > 1)
                    rc |= 8 if bad else 0
                    out.append(f"{k}/{len(names)} Test #{k}: {t} {'.' * 30} {'***Failed ' if bad else '  Passed'}"
                               f"    0.01 sec")
                return subprocess.CompletedProcess(cmd, rc, "\n".join(out) + "\n", "")

            buf = io.StringIO()
            with mock.patch.object(mutation_check, "MUTATIONS", mutations), \
                    mock.patch.object(mutation_check.subprocess, "run", run), \
                    mock.patch.object(sys, "argv", ["mutation_check.py", str(src), str(src / "build"), sel]), \
                    contextlib.redirect_stdout(buf):
                rc = mutation_check.main()
            return rc, buf.getvalue(), calls

    M = [("f.cpp", "x = 1", "x = 2", ["t_a"])]

    def test_kill(self):
        rc, out, calls = self.run_main(self.M, ["t_a"], lambda t, body, n, par: "x = 2" in body)
        self.assertEqual(rc, 0, out)
        self.assertIn("KILLED t_a", out)
        self.assertEqual(len(calls), 3)  # the unmutated source, the mutant, the restored source

    def test_failing_before_the_mutant_is_no_kill(self):
        rc, out, _ = self.run_main(self.M, ["t_a"], lambda t, body, n, par: True)
        self.assertEqual(rc, 1)
        self.assertIn("SURVIVED t_a", out)
        self.assertIn("t_a does not pass on the unmutated source", out)

    def test_restored_source_checked_again(self):
        # t_a passes on the unmutated source and fails under the mutant, then fails on the restored source too (a
        # stale object, a scratch file left behind): its kill is not trusted
        rc, out, _ = self.run_main(self.M, ["t_a"], lambda t, body, n, par: "x = 2" in body or n >= 3)
        self.assertEqual(rc, 1)
        self.assertIn("t_a does not pass on the restored source", out)

    def test_one_test_at_a_time(self):
        # t_a and t_b share a scratch file: run together they both fail (CTEST_PARALLEL_LEVEL = 8 is set above)
        m = [("f.cpp", "x = 1", "x = 2", ["t_a", "t_b"])]
        rc, out, calls = self.run_main(m, ["t_a", "t_b"], lambda t, body, n, par: par or "x = 2" in body)
        self.assertEqual(rc, 0, out)
        self.assertIn("KILLED t_a", out)
        self.assertIn("KILLED t_b", out)

    def test_names_matched_literally(self):
        m = [("f.cpp", "x = 1", "x = 2", ["t+a"])]
        rc, out, _ = self.run_main(m, ["t+a", "ttta"], lambda t, body, n, par: "x = 2" in body)
        self.assertEqual(rc, 0, out)
        self.assertIn("KILLED t+a", out)

    def test_bad_selection_runs_nothing(self):
        with self.assertRaises(SystemExit) as e:
            self.run_main(self.M * 3, ["t_a"], lambda t, body, n, par: False, sel="1,2,")
        self.assertIn("'1,2,'", str(e.exception.code))


if __name__ == "__main__":
    unittest.main()
