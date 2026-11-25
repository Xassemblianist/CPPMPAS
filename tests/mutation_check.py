#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Show that every test can fail: apply one deliberate one-line bug at a time to the source, rebuild, run the tests
it must break, require each of them to fail, restore the source. Run it through the workstation gate.

usage: mutation_check.py <source dir> <build dir> [first]   (first: index of the first mutant to apply, default 0)
"""
import re
import subprocess
import sys
from pathlib import Path

R, W = "src/io/mpas_reader.cpp", "src/io/mpas_writer.cpp"
A, V = "src/core_atmosphere/atm_init.cpp", "src/operators/vector_operations.cpp"
# (file, original line fragment, mutated fragment, tests that must fail)
MUTATIONS = [
    ("src/core/decomp.cpp", "    for (int j = 1; j <= nEdgesOnCell[i - 1]; ++j) {\n      const int e",
     "    for (int j = nEdgesOnCell[i - 1]; j >= 1; --j) {\n      const int e",
     ["unit_decomp", "order_MW", "order_SC", "order_JW", "order_init_MW"]),
    ("src/core/decomp.cpp", "        field(j, i) = sortedID[k - 1].second;", "        field(j, i) = sortedID[k - 1].first;",
     ["unit_decomp", "order_MW", "order_SC", "order_JW", "order_init_MW"]),
    (W, "      for (int j = m.nEdgesOnCell(i) + 1; j <= m.maxEdges; ++j) put(j, i, 0);",
     "      for (int j = m.nEdgesOnCell(i) + 1; j <= m.maxEdges; ++j) put(j, i, id(a(j, i)));", ["roundtrip_grid_JW"]),
    (W, "    local_to_file(o, inner, local_values.data(), out.data());",
     "    file_to_local(o, inner, local_values.data(), out.data());",
     ["roundtrip_grid_MW", "roundtrip_grid_SC", "roundtrip_grid_JW", "roundtrip_init_MW", "roundtrip_init_SC",
      "history_MW", "history_SC"]),
    (W, "  NcFile f = NcFile::create(path, NC_64BIT_DATA);", "  NcFile f = NcFile::create(path, NC_64BIT_OFFSET);",
     ["history_MW", "history_SC"]),
    (W, "      f.put_att_text(n, k, v);", "      f.put_att_text(n, k == \"long_name\" ? \"description\" : k, v);",
     ["history_MW", "history_SC"]),
    (R, "  if (m.is_periodic)\n", "  if (false && m.is_periodic)\n", ["err_periodic"]),
    (R, "    if (!f.has_var(n)) throw Error(f.path() + \": required mesh variable \"",
     "    if (false) throw Error(f.path() + \": required mesh variable \"", ["err_missing_var"]),
    (R, "  if (vi.type != want) {", "  if (false) {", ["err_float_var"]),
    (R, "fmt != NC_FORMAT_CDF5)", "fmt != NC_FORMAT_CDF5 && fmt != NC_FORMAT_NETCDF4)", ["err_netcdf4"]),
    (R, "      if (v[i] != static_cast<int>(i) + 1)", "      if (false)", ["err_index_ids"]),
    (R, "      if (cellsOnCell[(j - 1) + m.maxEdges * (i - 1)] == 0)",
     "      if (cellsOnCell[(j - 1) + m.maxEdges * (i - 1)] == -1)", ["err_boundary"]),
    (R, "      if (a->data()[k] != 0)", "      if (false)", ["err_bdymask"]),
    (R, "      if (vi.shape[k] != 1)", "      if (false)", ["err_two_records"]),
    (R, "    } else if (init) {", "    } else if (false) {", ["err_missing_init_var"]),
    (R, "    if (!f.has_var(n)) throw Error(f.path() + \": active scalar \"",
     "    if (false) throw Error(f.path() + \": active scalar \"", ["err_missing_scalar"]),
    (R, "  if (vi.dims != want_dims)", "  if (false)", ["err_bad_dims"]),
    (R, "    if (std::find(order.begin(), order.end(), s) == order.end())", "    if (false)", ["err_unknown_scalar"]),
    (R, "  throw Error(\"variable '\" + name + \"' has no CPPMPAS port (not read or written)\");",
     "  throw Error(\"variable '\" + name + \"' was not read\");", ["err_unported_field"]),
    ("src/io/registry.cpp", "  if (!r) throw Error(\"variable '\" + name + \"' is not defined in MPAS-Atmosphere's Registry\");",
     "  if (!r) throw Error(\"variable '\" + name + \"' is unknown\");", ["err_unknown_field"]),
    # --- OPM03: reader on the reference files, atm_mpas_init_block against the references, libm ports
    (R, "    if (kInputNotUsed.count(n)) continue;", "    if (false) continue;",
     ["readrestart_MW", "readrestart_SC", "readrestart_JW", "roundtrip_refinit_MW", "roundtrip_refinit_SC",
      "roundtrip_refinit_JW", "order_refinit_MW", "order_refinit_SC", "order_refinit_JW"]),
    (W, "    local_to_file(o, inner, local_values.data(), out.data());",
     "    file_to_local(o, inner, local_values.data(), out.data());",
     ["roundtrip_refinit_MW", "roundtrip_refinit_SC", "roundtrip_refinit_JW", "init_restart_MW"]),  # cells keep file order
    ("src/core/libm.hpp", "    *sinx = g::copysign(g::do_cos(a, da), x);", "    *sinx = g::do_cos(a, da);",
     ["unit_libm"]),
    ("src/core/libm.hpp", "inline double mvec_pow(double x, double y) { return ::cppwrf::glibc::mvec::pow(x, y); }",
     "inline double mvec_pow(double x, double y) { return ::cppwrf::glibc::pow(x, y); }", ["unit_libm"]),
    (A, "      d.exner(k, iCell) = vector ? libm::mvec_pow(xe, rcv) : libm::pow(xe, rcv);",
     "      d.exner(k, iCell) = libm::pow(xe, rcv);", ["init_step0_MW", "init_step0_SC", "init_step0_JW"]),
    (A, "      d.exner_base(k, iCell) = vector ? libm::mvec_pow(xb, rcv) : libm::pow(xb, rcv);",
     "      d.exner_base(k, iCell) = libm::pow(xb, rcv);", ["init_restart_MW", "init_restart_SC", "init_restart_JW"]),
    (A, "      rho_zz(k, iCell) = s.rho(k, iCell) / zz(k, iCell);", "      rho_zz(k, iCell) = s.rho(k, iCell) * (1.0 / zz(k, iCell));",
     ["init_step0_MW"]),  # SC: zz = 1
    (A, "      theta_m(k, iCell) = s.theta(k, iCell) * (1. + rvord * s.scalars(index_qv, k, iCell));",
     "      theta_m(k, iCell) = s.theta(k, iCell) + s.theta(k, iCell) * rvord * s.scalars(index_qv, k, iCell);",
     ["init_step0_SC"]),
    (A, "    for (int k = 1; k <= nVertLevels; ++k) d.rtheta_base(k, iCell) = s.theta_base(k, iCell) * s.rho_base(k, iCell);",
     "    for (int k = 1; k <= nVertLevels; ++k) d.rtheta_base(k, iCell) = s.rho_base(k, iCell) * s.theta_base(k, iCell) * 1.0000000000000002;",
     ["init_restart_MW"]),
    (A, "        m.dss(k, iCell) = config_xnutr * (s * s);", "        m.dss(k, iCell) = config_xnutr * s * s;",
     ["init_restart_MW"]),
    (A, "        const real s = libm::sin(0.5 * pii * (z - config_zd) / (zt - config_zd));",
     "        const real s = libm::sin(0.5 * pii * ((z - config_zd) / (zt - config_zd)));", ["init_restart_MW"]),
    (V, "  xin[2] = xin[2] / mag;", "  xin[2] = xin[2] * (1.0 / mag);",
     ["init_restart_JW"]),  # MW, SC planar: xin[2] = 0
    (V, "      mpas_unit_vec_in_r3(v);\n      for (int d = 1; d <= 3; ++d) lvu(d, iCell) = v[d - 1];",
     "      for (int d = 1; d <= 3; ++d) lvu(d, iCell) = v[d - 1];", ["init_restart_JW"]),
    (V, "    for (int d = 1; d <= 3; ++d) normalDotRHat = normalDotRHat + env(d, iEdge) * rHat[d - 1];",
     "    for (int d = 3; d >= 1; --d) normalDotRHat = normalDotRHat + env(d, iEdge) * rHat[d - 1];",
     ["init_restart_JW"]),  # MW, SC planar: two of the three terms are 0
    (V, "    alpha = alpha / static_cast<real>(pointCount);", "    alpha = alpha * (1.0 / static_cast<real>(pointCount));",
     ["init_restart_MW", "init_restart_JW"]),
    (V, "      for (int d = 1; d <= 3; ++d) s = s + SP(i, d) * PB(r, d);",
     "      for (int d = 3; d >= 1; --d) s = s + SP(i, d) * PB(r, d);", ["init_restart_JW"]),
    (V, "    X[I - 1] = X[I - 1] / a(INDX[I - 1], I);", "    X[I - 1] = X[I - 1] * (1.0 / a(INDX[I - 1], I));",
     ["init_restart_MW"]),
    ("src/core/namelist.cpp", "  if (c.config_do_restart)\n", "  if (false)\n", ["err_do_restart"]),
    ("src/core/namelist.cpp", "    if (t == \"/\") { in_group = false; continue; }", "    if (false) { in_group = false; continue; }",
     ["init_run_MW", "init_run_SC", "init_run_JW"]),
]


def main():
    src, build = Path(sys.argv[1]), Path(sys.argv[2])
    first = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    covered, bad = set(), []
    for path, old, new, tests in MUTATIONS[first:]:
        p = src / path
        orig = p.read_text()
        if orig.count(old) != 1:
            bad.append(f"{path}: fragment found {orig.count(old)} times: {old!r}")
            continue
        try:
            p.write_text(orig.replace(old, new))
            b = subprocess.run(["cmake", "--build", str(build), "-j4"], capture_output=True, text=True)
            if b.returncode != 0:
                bad.append(f"{path}: mutant does not build")
                continue
            for t in tests:
                r = subprocess.run(["ctest", "--test-dir", str(build), "-R", f"^{t}$", "-FA", "inputs"],
                                   capture_output=True, text=True)
                # the named test's own line in the summary (its fixture setup, if any, runs first)
                ran = re.search(rf"Test +#\d+: {re.escape(t)} \.+", r.stdout) is not None
                failed = r.returncode != 0 and re.search(rf"^\s*\d+ - {re.escape(t)} \(Failed\)$", r.stdout, re.M) is not None
                print(f"{'KILLED' if failed else 'SURVIVED'} {t:22s} by {path}: {new.strip()[:70]}")
                if failed:
                    covered.add(t)
                else:
                    bad.append(f"{t} survived {path}: {new.strip()} (ran={ran})")
        finally:
            p.write_text(orig)
    subprocess.run(["cmake", "--build", str(build), "-j4"], capture_output=True, text=True, check=True)
    for b in bad:
        print("PROBLEM", b)
    print(f"{len(covered)} tests killed by at least one mutant")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
