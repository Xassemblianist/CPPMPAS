#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""tests/mutation_check.py's reading of CTest output: a test that failed or crashed kills the mutant; a test that
passed, or was not run because its fixture setup failed (CTest's "***Not Run"), does not.

usage: test_mutation_check.py
"""
import sys
import unittest
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import mutation_check  # noqa: E402

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

    def test_absent_or_prefix(self):
        self.assertFalse(self.k(FAILED, "initdump_JW"))  # a prefix of a test name is not that test
        self.assertFalse(self.k("", "initdump_JW_ke"))


if __name__ == "__main__":
    unittest.main()
