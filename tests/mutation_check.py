#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Show that every test can fail: apply one deliberate one-line bug at a time to the source, rebuild, run the tests
it must break, require each of them to fail, restore the source. Run it through the workstation gate.

usage: mutation_check.py <source dir> <build dir> [first | i,j,...]   (index of the first mutant to apply, default 0,
       or a comma-separated list of mutant indices)
"""
import re
import subprocess
import sys
from pathlib import Path

R, W = "src/io/mpas_reader.cpp", "src/io/mpas_writer.cpp"
N, X = "src/io/netcdf_file.cpp", "tools/cppmpas.cpp"
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
    (R, "      if (v == 0 && boundary_element)", "      if (false && boundary_element)",
     ["err_boundary", "err_boundary_edge", "err_boundary_vertex"]),
    (R, "      if (a->data()[k] != 0)", "      if (false)", ["err_bdymask"]),
    (R, "      if (vi.shape[k] != 1)", "      if (false)", ["err_two_records"]),
    (R, "    } else if (init) {", "    } else if (false) {", ["err_missing_init_var"]),
    (R, "    if (!f.has_var(n)) throw Error(f.path() + \": active scalar \"",
     "    if (false) throw Error(f.path() + \": active scalar \"", ["err_missing_scalar"]),
    (R, "  if (vi.dims != want_dims)", "  if (false)", ["err_bad_dims"]),
    (R, "    if (std::find(order.begin(), order.end(), *it) == order.end())", "    if (false)", ["err_unknown_scalar"]),
    (R, "  throw Error(\"variable '\" + name + \"' has no CPPMPAS port (not read or written)\");", "  return r;",
     ["err_unported_field", "unit_io_core"]),
    # a name the Registry does not define: write_history meets registry_var's check first and field_ref's second, so
    # each guard is tested on its own (unit_io_core); err_unknown_field checks the CLI message
    ("src/io/registry.cpp", "  if (!r) throw Error(\"variable '\" + name + \"' is not defined in MPAS-Atmosphere's Registry\");",
     "  if (false) throw Error(\"variable '\" + name + \"' is not defined in MPAS-Atmosphere's Registry\");",
     ["unit_io_core"]),
    (R, "  if (!find_registry_var(name)) throw", "  if (false) throw", ["unit_io_core"]),
    # --- reader on the reference files, atm_mpas_init_block against the references, libm ports
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
    # --- review fixes: index arithmetic, NetCDF wrapper, reader checks, writer layout, CLI, build flags
    ("src/core/field.hpp", "    return D(i - 1) + D(e0) * (D(j - 1) + D(e1) * (D(k - 1) + D(e2) * D(l - 1)));",
     "    return (i - 1) + e0 * ((j - 1) + e1 * ((k - 1) + e2 * (l - 1)));", ["unit_io_core"]),
    ("src/core/mesh.hpp", "  int n() const { return local_to_file.empty() ? 0 : static_cast<int>(local_to_file.size()) - 1; }",
     "  int n() const { return static_cast<int>(local_to_file.size()) - 1; }", ["unit_io_core"]),
    (N, "long record) {\n  std::vector<std::size_t> s, c;\n  record_range(name, record, s, c);  // record >= 0 needs a variable whose first dimension is the record dimension\n  if (product(c) != v.size()) throw Error(path_ + \": write \" + name + \": size mismatch\");\n  check(nc_put_vara_double(",
     "long record) {\n  std::vector<std::size_t> s, c;\n  record_range(name, -1, s, c); if (record >= 0) { s[0] = record; c[0] = 1; }\n  if (product(c) != v.size()) throw Error(path_ + \": write \" + name + \": size mismatch\");\n  check(nc_put_vara_double(",
     ["unit_io_core"]),
    (N, "  if (st == NC_EBADDIM) return false;", "  if (st != NC_NOERR) return false;", ["unit_io_core"]),
    (N, "  f.check(nc_set_fill(id, NC_NOFILL, &old), \"set_fill\");", "  (void)old;", ["unit_io_core"]),
    (R, "    check_count(f, \"nEdgesOnEdge\", nEdgesOnEdge, 0, m.maxEdges2, \"maxEdges2\");", "", ["err_nEdgesOnEdge"]),
    (R, "    check_count(f, \"nEdgesOnCell\", nEdgesOnCell, 1, m.maxEdges, \"maxEdges\");", "", ["err_nEdgesOnCell"]),
    (R, "      if (v < 1 || v > n)", "      if (false)",
     ["err_cellsOnEdge_range", "err_cellsOnCell_neg", "err_verticesOnEdge_zero", "err_edgesOnVertex_range",
      "err_edgesOnEdge_range"]),
    (R, "    const int nu = used ? (*used)[i - 1] : inner;", "    const int nu = inner;", ["order_JW", "roundtrip_grid_JW"]),
    (R, "  if (n < 1 || n > max)", "  if (n > max)", ["err_zero_dim"]),
    (R, "    for (const auto& d : rv.dims) vertical |= (d == \"nVertLevels\" || d == \"nVertLevelsP1\");",
     "    for (const auto& d : rv.dims) vertical |= (d == \"nVertLevels\");", ["grid_read_of_init"]),
    (R, "    if (std::find(scalars.begin(), it, *it) != it) throw", "    if (false) throw", ["err_dup_scalar"]),
    (R, "  if (std::find(scalars.begin(), scalars.end(), \"qv\") == scalars.end())", "  if (false)", ["err_no_qv"]),
    (R, "(rtrim(file_id) + \"\\n\" + rtrim(parent_id)).substr(0, kStrKIND);", "(rtrim(file_id) + \"\\n\" + rtrim(parent_id));",
     ["global_atts"]),
    (R, "    in.warnings.push_back(\"Attribute file_id not found;", "    if (false) in.warnings.push_back(\"Attribute file_id not found;",
     ["global_atts"]),
    (R, "  if (a.len > kStrKIND)", "  if (false)", ["err_long_att"]),
    (W, "  for (auto it = dims.rbegin(); it != dims.rend(); ++it) {", "  for (auto it = dims.begin(); it != dims.end(); ++it) {",
     ["history_MW", "history_SC", "history_JW"]),
    (W, "opt.file_id.empty() ? random_file_id() : opt.file_id", "opt.file_id", ["history_SC", "history_JW"]),
    (W, "    if (ref.c->empty()) throw", "    if (false) throw", ["unit_io_core"]),
    (W, "  for (const auto& a : in.source.global_atts) f.put_att(\"\", a);",
     "  for (const auto& a : in.source.global_atts) if (a.name != \"parent_id\") f.put_att(\"\", a);",
     ["roundtrip_init_MW", "roundtrip_init_SC", "roundtrip_refinit_MW", "roundtrip_refinit_SC", "roundtrip_refinit_JW"]),
    (X, "    if (r.scalar_slot > 0) continue;", "    if (r.c || r.scalar_slot > 0) continue;",
     ["order_init_MW", "order_init_SC", "order_refinit_MW", "order_refinit_SC", "order_refinit_JW"]),
    (X, "f.put_double(\"scalars\", std::vector<double>(sc.data(), sc.data() + sc.size()));",
     "f.put_double(\"scalars\", std::vector<double>(sc.size(), 0.0));",
     ["order_init_MW", "order_init_SC", "order_refinit_SC"]),  # MW, JW reference init files: qv = 0
    (X, "    if (x.empty()) throw", "    if (x.empty()) continue; else if (false) throw", ["err_empty_name"]),
    (X, "  if (v.empty() || s.back() == ',')", "  if (false)", ["err_empty_list", "err_trailing_comma"]),
    ("CMakeLists.txt", "add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-ffp-contract=off>)\n",
     "add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-ffp-contract=off>)\nset(CMAKE_CXX_FLAGS_RELEASE \"-O3\")\n",
     ["cmake_rejects_Ofast", "cmake_rejects_ffast-math"]),
    ("tools/gen_registry_table.py", "        return \"\"\n    if mv == \"FILLVAL\"", "        return \"\"\n    return \"\"\n    if mv == \"FILLVAL\"",
     ["registry_table_current"]),
    ("tools/mpas_compare.py", "if ra.tobytes() != rb.tobytes():", "if False:", ["compare_unit"]),
    ("src/core/decomp.cpp", "    for (int j = 1; j <= nEdgesOnCell[i - 1]; ++j) h.insert(",
     "    for (int j = 1; j <= maxEdges; ++j) h.insert(", ["unit_decomp"]),
    ("src/core/decomp.cpp", "  int result = n2 + 1;", "  int result = n2;", ["unit_decomp"]),
]


def main():
    src, build = Path(sys.argv[1]), Path(sys.argv[2])
    sel = sys.argv[3] if len(sys.argv) > 3 else "0"
    chosen = [MUTATIONS[int(i)] for i in sel.split(",")] if "," in sel else MUTATIONS[int(sel):]
    covered, bad = set(), []
    for path, old, new, tests in chosen:
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
                # "(Failed)", or "(Subprocess aborted)", "(SEGFAULT)", ... for a test the mutant makes crash
                failed = r.returncode != 0 and re.search(rf"^\s*\d+ - {re.escape(t)} \(.+\)$", r.stdout, re.M) is not None
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
