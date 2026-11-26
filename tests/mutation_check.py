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
NL, L = "src/core/namelist.cpp", "src/core/libm.hpp"
S, RC = "src/core_atmosphere/atm_solve_diagnostics.cpp", "src/operators/vector_reconstruct.cpp"


def D(fields, cases="MW SC JW"):
    """initdump tests (tests/CMakeLists.txt) of the given fields and cases"""
    return [f"initdump_{c}_{f}" for c in cases.split() for f in fields.split()]


ALL3 = "MW SC JW"
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
    (W, "    for (const auto& [k, v] : rv.attrs) f.put_att_text(n, k, v);",
     "    for (const auto& [k, v] : rv.attrs) f.put_att_text(n, k == \"long_name\" ? \"description\" : k, v);",
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
    ("src/core/libm.hpp", "  return ::cppwrf::glibc::mvec::pow(arg(x), arg(y));",
     "  return ::cppwrf::glibc::pow(arg(x), arg(y));", ["unit_libm"]),
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
    (NL, "    if (!ok) {", "    if (false) {", ["namelist_MW"]),
    (NL, "    if (c == '/') {", "    if (c == '#') {", ["namelist_MW", "init_run_MW", "init_run_SC", "init_run_JW"]),
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
    (W, "    if (!held) throw", "    if (false) throw", ["unit_io_core", "err_unheld_field"]),
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
    # compare_mutants fails when a unit test of the comparator is weakened so that one of its mutants survives
    ("tests/compare/test_mpas_compare.py", "        self.assertEqual(code, 2)", "        self.assertIn(code, (0, 1, 2))",
     ["compare_mutants"]),
    ("src/core/decomp.cpp", "    for (int j = 1; j <= nEdgesOnCell[i - 1]; ++j) h.insert(",
     "    for (int j = 1; j <= maxEdges; ++j) h.insert(", ["unit_decomp"]),
    ("src/core/decomp.cpp", "  int result = n2 + 1;", "  int result = n2;", ["unit_decomp"]),
    # --- review fixes of the init block: namelist, atm_core_init, guards, unverified switch, writer, tests, libm
    (NL, "    if (!find_option(record, name))", "    if (false)", ["namelist_MW"]),
    (NL, "    if (kv.count(name)) throw", "    if (false) throw", ["namelist_MW"]),
    (NL, "    if (c == '&' || c == '$') throw", "    if (false) throw", ["namelist_MW"]),
    (NL, "      if (k == j || k >= s.size() || s[k] != '=')", "      if (false)", ["namelist_MW"]),
    (NL, "      if (v.text.find('*') != std::string::npos) throw", "      if (false) throw", ["namelist_MW"]),
    (NL, "    if (!std::regex_match(v.text, kReal)) throw", "    if (false) throw", ["namelist_MW"]),
    (NL, "    if (!std::regex_match(v.text, kInt)) throw", "    if (false) throw", ["namelist_MW"]),
    (NL, "    if (!v.quoted) throw", "    if (false) throw", ["namelist_MW"]),
    (NL, "  if (v.quoted) throw", "  if (false) throw", ["namelist_MW"]),
    (NL, "      if (!is_record(record)) throw", "      if (false) throw", ["namelist_MW"]),
    (NL, "      if (c != '&') throw", "      if (false) throw", ["namelist_MW"]),
    (NL, "  if (!record.empty()) throw", "  if (false) throw", ["namelist_MW"]),
    (NL, "  if (l == \"f\" || l == \".f.\" ||", "  if (l == \"f\" ||", ["namelist_MW"]),
    (NL, "  r(\"config_len_disp\", c.config_len_disp);", "", ["init_run_MW", "init_run_SC", "init_run_JW"]),
    (NL, "  t(\"config_start_time\", c.config_start_time);", "", ["init_step0_MW", "init_step0_SC", "init_step0_JW"]),
    ("tools/gen_namelist_table.py", "            if v not in acc:", "            if True:", ["namelist_table_current"]),
    # the reference init files already hold xtime = config_start_time: the init_step0_* comparison cannot see it
    (A, "  s.xtime = stamp;", "", ["unit_init_guards"]),
    (A, "    nominalMinDc = cfg.config_len_disp;", "", ["unit_init_guards", "init_restart_MW", "init_restart_SC"]),
    (A, "    cfg.config_len_disp = nominalMinDc;", "", ["unit_init_guards"]),
    (A, "  if (!std::regex_match(cfg.config_start_time, f, canonical))", "  if (false)", ["unit_init_guards"]),
    (A, "  if (!cfg.config_apply_lbcs && maxmask > 0)", "  if (false)", ["unit_init_guards"]),
    (A, "  if (cfg.config_apply_lbcs && maxmask == 0)", "  if (false)", ["unit_init_guards"]),
    (A, "  if (cfg.config_epssm != 0.0)", "  if (false)", ["unit_init_guards"]),
    (A, "  if (cfg.config_positive_definite)", "  if (false)", ["unit_init_guards"]),
    (A, "    if (n > 25) throw", "    if (n > 26) throw", ["unit_init_guards"]),
    (A, "  if (m.is_periodic)\n    throw Error(\"atm_initialize_deformation_weights", "  if (false)\n    throw Error(\"atm_initialize_deformation_weights",
     ["unit_init_guards"]),
    (A, "      if (a->data()[k] != 0) throw Error(\"mpas_atm_setup_bdy_masks", "      if (false) throw Error(\"mpas_atm_setup_bdy_masks",
     ["unit_init_guards"]),
    (V, "    if (m.cellsOnEdge(1, iEdge) == nCells + 1 || m.cellsOnEdge(2, iEdge) == nCells + 1)",
     "    if (m.cellsOnEdge(1, iEdge) == nCells + 1)", ["unit_init_guards"]),
    (R, "    buf = f.get_double(rv.name, record);", "    buf = f.get_double(rv.name, record), buf[0] = buf[0] * 0.5;",
     ["readrestart_MW", "readrestart_SC", "readrestart_JW"]),
    ("tests/compare_arrays.py", "            out.append((k, arr.dtype.str, arr.tobytes()))", "            out.append((k, \"x\", str(a)))",
     ["compare_arrays_unit"]),
    ("tests/compare_arrays.py", "            if attributes(x) != attributes(y):", "            if False:", ["compare_arrays_unit"]),
    (L, "double cos(double x) { return ::cppwrf::glibc::cos(arg(x)); }",
     "double cos(double x) { return ::cppwrf::glibc::mvec::cos(arg(x)); }", ["unit_libm"]),
    (L, "    *sinx = g::copysign(g::do_cos(a, da), x);", "    *sinx = g::copysign(g::do_cos(y, hp1), x);", ["unit_libm"]),
    (L, "1.0 ? __builtin_nan(\"\") : ::cppwrf::glibc::asin(x);", "1.0 ? -__builtin_nan(\"\") : ::cppwrf::glibc::asin(x);",
     ["unit_libm"]),
    # --- OPM05: every initdump field (atm_mpas_init_block at t = 0) against mpas842-dbl-mw-003 / sc-002 / jw-002
    # atm_compute_solve_diagnostics
    (S, "h_edge(k, iEdge) = 0.5 * (h(k, cell1) + h(k, cell2));", "h_edge(k, iEdge) = 0.5 * (h(k, cell1) + h(k, cell1));",
     D("rho_edge", "MW SC JW")),
    (S, "ke_edge(k, iEdge) = efac * (u(k, iEdge) * u(k, iEdge));", "ke_edge(k, iEdge) = (efac * u(k, iEdge)) * u(k, iEdge);",
     D("ke", "SC JW")),
    (S, "    for (int i = 1; i <= vertexDegree; ++i) {", "    for (int i = vertexDegree; i >= 1; --i) {",
     D("vorticity pv_vertex", "SC JW")),
    (S, "divergence(k, iCell) = divergence(k, iCell) * r;", "divergence(k, iCell) = divergence(k, iCell) / m.areaCell(iCell);",
     D("divergence", "SC JW")),
    (S, "(ke_edge(k, e1) + ke_edge(k, e2) + ke_edge(k, e3)) * r;", "(ke_edge(k, e3) + ke_edge(k, e2) + ke_edge(k, e1)) * r;",
     D("ke", "SC JW")),
    (S, "(1. - ke_fact) * m.kiteAreasOnVertex(j, iVertex) * ke_vertex(k, iVertex) * r;",
     "(1. - ke_fact) * m.kiteAreasOnVertex(j, iVertex) * (ke_vertex(k, iVertex) * r);", D("ke", "SC JW")),
    (S, "      for (int i = 1; i <= m.nEdgesOnEdge(iEdge); ++i) {", "      for (int i = m.nEdgesOnEdge(iEdge); i >= 1; --i) {",
     D("v", "MW SC JW")),
    (S, "pv_edge(k, iEdge) = 0.5 * (pv_vertex(k, m.verticesOnEdge(1, iEdge)) + pv_vertex(k, m.verticesOnEdge(2, iEdge)));",
     "pv_edge(k, iEdge) = 0.5 * (pv_vertex(k, m.verticesOnEdge(1, iEdge)) + pv_vertex(k, m.verticesOnEdge(1, iEdge)));",
     D("pv_edge", "MW SC JW")),
    (S, "pv_cell(k, iCell) + m.kiteAreasOnVertex(j, iVertex) * pv_vertex(k, iVertex) * r;",
     "pv_cell(k, iCell) + m.kiteAreasOnVertex(j, iVertex) * (pv_vertex(k, iVertex) * r);", D("pv_cell gradPVn", "SC JW")),
    (S, "      gradPVt(k, iEdge) = (pv_vertex(k, v2) - pv_vertex(k, v1)) * r1;",
     "      gradPVt(k, iEdge) = (pv_vertex(k, v2) - pv_vertex(k, v1)) * r2;", D("gradPVt", "MW SC JW")),
    (S, "  const real r = config_apvm_upwinding * dt;", "  const real r = config_apvm_upwinding;", D("pv_edge", "MW SC JW")),
    (S, "  if (!(config_apvm_upwinding > 0.0))", "  if (!(config_apvm_upwinding > 1.0))",
     ["initdump_run_MW", "initdump_run_SC", "initdump_run_JW"]),
    # MW: u is horizontally uniform (vorticity = +-1.26e-17, pv_cell = gradPVn = 0, ke = 50, divergence = 1.6e-18),
    # so the reorderings above are exact there; these bugs change the MW values
    (S, "      const real s = m.edgesOnVertex_sign(i, iVertex) * m.dcEdge(iEdge);",
     "      const real s = m.edgesOnVertex_sign(i, iVertex) * m.dvEdge(iEdge);", D("vorticity pv_vertex", "MW SC JW")),
    (S, "      const real s = m.edgesOnCell_sign(i, iCell) * m.dvEdge(iEdge);",
     "      const real s = m.edgesOnCell_sign(i, iCell) * m.dcEdge(iEdge);", D("divergence", "MW SC JW")),
    (S, "ke(k, iCell) = ke(k, iCell) + 0.25 * ke_edge(k, iEdge);", "ke(k, iCell) = ke(k, iCell) + 0.5 * ke_edge(k, iEdge);",
     D("ke", "MW SC JW")),
    (S, "pv_cell(k, iCell) = 0.0;\n    const real r = m.invAreaCell(iCell);\n    for (int i = 1; i <= m.nEdgesOnCell(iCell); ++i) {\n      const int iVertex = m.verticesOnCell(i, iCell);",
     "pv_cell(k, iCell) = 0.0;\n    const real r = m.invAreaCell(iCell);\n    for (int i = 1; i <= m.nEdgesOnCell(iCell); ++i) {\n      const int iVertex = m.verticesOnCell(1, iCell);",
     D("pv_cell gradPVn", "SC JW") + D("pv_cell", "MW SC")),  # MW gradPVn: neighbouring cells change alike
    (S, "      gradPVn(k, iEdge) = (pv_cell(k, c2) - pv_cell(k, c1)) * r2;",
     "      gradPVn(k, iEdge) = (pv_vertex(k, c2) - pv_vertex(k, c1)) * r2;", D("gradPVn", "MW SC JW")),
    (S, "  if (!(config_apvm_upwinding > 0.0))", "  if (false)", ["unit_init_guards"]),
    (A, "    if (nVertLevels % 2 != 0 && std::find(pow_tail.begin(), pow_tail.end(), n) != pow_tail.end()) continue;",
     "    if (false && std::find(pow_tail.begin(), pow_tail.end(), n) != pow_tail.end()) continue;", ["unit_init_guards"]),
    # mpas_reconstruct
    (RC, "uReconstructX(k, iCell) + coeffs_reconstruct(1, i, iCell) * u(k, iEdge);",
     "uReconstructX(k, iCell) + coeffs_reconstruct(2, i, iCell) * u(k, iEdge);", D("uReconstructX uReconstructZonal", "MW SC JW")),
    (RC, "uReconstructY(k, iCell) + coeffs_reconstruct(2, i, iCell) * u(k, iEdge);",
     "uReconstructY(k, iCell) + coeffs_reconstruct(1, i, iCell) * u(k, iEdge);", D("uReconstructY uReconstructMeridional", "MW SC JW")),
    (RC, "uReconstructZ(k, iCell) + coeffs_reconstruct(3, i, iCell) * u(k, iEdge);",
     "uReconstructZ(k, iCell) + coeffs_reconstruct(1, i, iCell) * u(k, iEdge);", D("uReconstructZ", "MW SC JW")),
    (RC, "        uReconstructZonal(k, iCell) = uReconstructX(k, iCell);", "        uReconstructZonal(k, iCell) = uReconstructY(k, iCell);",
     D("uReconstructZonal", "MW SC")),
    (RC, "        uReconstructMeridional(k, iCell) = uReconstructY(k, iCell);",
     "        uReconstructMeridional(k, iCell) = uReconstructX(k, iCell);", D("uReconstructMeridional", "MW SC")),
    # sincos versus sin / cos of latCell (VR:310-311): no difference on any JW cell, an equivalent mutant there
    (RC, "            -(uReconstructX(k, iCell) * clon + uReconstructY(k, iCell) * slon) * slat + uReconstructZ(k, iCell) * clat;",
     "            -(uReconstructX(k, iCell) * clon * slat + uReconstructY(k, iCell) * slon * slat) + uReconstructZ(k, iCell) * clat;",
     D("uReconstructMeridional", "JW")),
    (RC, "      libm::sincos(m.lonCell(iCell), &slon, &clon);", "      libm::sincos(m.lonCell(iCell), &clon, &slon);",
     D("uReconstructZonal uReconstructMeridional", "JW")),
    # the arrays OPM03 could not compare
    (A, "m.invAreaCell(iCell) = 1.0 / m.areaCell(iCell);",
     "m.invAreaCell(iCell) = 1.0 / std::sqrt(m.areaCell(iCell)) / std::sqrt(m.areaCell(iCell));", D("invAreaCell", "MW SC JW")),
    (A, "m.invDvEdge(iEdge) = 1.0 / m.dvEdge(iEdge);",
     "m.invDvEdge(iEdge) = 1.0 / std::sqrt(m.dvEdge(iEdge)) / std::sqrt(m.dvEdge(iEdge));", D("invDvEdge", "SC JW")),
    (A, "m.invDcEdge(iEdge) = 1.0 / m.dcEdge(iEdge);",
     "m.invDcEdge(iEdge) = 1.0 / std::sqrt(m.dcEdge(iEdge)) / std::sqrt(m.dcEdge(iEdge));", D("invDcEdge", "SC JW")),
    (A, "m.invAreaTriangle(iVertex) = 1.0 / m.areaTriangle(iVertex);",
     "m.invAreaTriangle(iVertex) = 1.0 / std::sqrt(m.areaTriangle(iVertex)) / std::sqrt(m.areaTriangle(iVertex));",
     D("invAreaTriangle", "MW SC JW")),
    (A, "m.invDvEdge(iEdge) = 1.0 / m.dvEdge(iEdge);", "m.invDvEdge(iEdge) = 1.0 / m.dcEdge(iEdge);", D("invDvEdge", "MW SC JW")),
    (A, "m.invDcEdge(iEdge) = 1.0 / m.dcEdge(iEdge);", "m.invDcEdge(iEdge) = 1.0 / m.dvEdge(iEdge);", D("invDcEdge", "MW SC JW")),
    (A, "          m.zb3_cell(k, i, iCell) = m.zb3(k, side, e);", "          m.zb3_cell(k, i, iCell) = m.zb(k, side, e);",
     D("zb3_cell", "MW SC JW")),
    (A, "        d.rw(k, iCell) = d.rw(k, iCell) - m.edgesOnCell_sign(i, iCell) *",
     "        d.rw(k, iCell) = d.rw(k, iCell) - ", D("rw", "MW SC JW")),
    (A, "        m.edgesOnCell_sign(i, iCell) = side == 1 ? 1.0 : -1.0;", "        m.edgesOnCell_sign(i, iCell) = 1.0;",
     D("edgesOnCell_sign", "MW SC JW")),
    (A, "          m.edgesOnVertex_sign(i, iVtx) = -1.0;", "          m.edgesOnVertex_sign(i, iVtx) = 1.0;",
     D("edgesOnVertex_sign", "MW SC JW")),
    (A, "          m.zb_cell(k, i, iCell) = m.zb(k, side, e);", "          m.zb_cell(k, i, iCell) = m.zb(k, 1, e);",
     D("zb_cell", "MW SC JW")),
    (A, "          m.zb3_cell(k, i, iCell) = m.zb3(k, side, e);", "          m.zb3_cell(k, i, iCell) = m.zb3(k, 1, e);",
     D("zb3_cell", "SC JW")),
    (A, "            m.kiteForCell(i, iCell) = j;", "            m.kiteForCell(i, iCell) = 1;", D("kiteForCell", "MW SC JW")),
    (A, "    for (int iCell = 1; iCell <= n; ++iCell) m.advCellsForEdge(iCell, iEdge) = cell_list[iCell];",
     "    for (int iCell = 1; iCell <= n; ++iCell) m.advCellsForEdge(iCell, iEdge) = cell_list[1];", D("advCellsForEdge", "MW SC JW")),
    (A, "    m.nAdvCellsForEdge(iEdge) = n;", "    m.nAdvCellsForEdge(iEdge) = n - 1;", D("nAdvCellsForEdge", "MW SC JW")),
    (A, "    j_in = slot(cell1);  // :1513-1517\n    adv(j_in, iEdge) = adv(j_in, iEdge) + 0.5;",
     "    j_in = slot(cell1);  // :1513-1517\n    adv(j_in, iEdge) = adv(j_in, iEdge) + 0.25;", D("adv_coefs", "MW SC JW")),
    (A, "m.adv_coefs_3rd.data()[k] = c * m.adv_coefs_3rd.data()[k];", "m.adv_coefs_3rd.data()[k] = m.adv_coefs_3rd.data()[k];",
     D("adv_coefs_3rd", "MW SC")),  # JW: config_coef_3rd_order = 1
    (A, "    adv3(j_in, iEdge) = adv3(j_in, iEdge) - m.deriv_two(1, 2, iEdge);",
     "    adv3(j_in, iEdge) = adv3(j_in, iEdge) + m.deriv_two(1, 2, iEdge);", D("adv_coefs_3rd", "MW SC JW")),
    # JW: config_h_ScaleWithMesh = true overwrites the fills; meshDensity = 1, so the bug must change the argument
    (A, "m.meshScalingDel2(iEdge) = 1.0 / libm::pow((m.meshDensity(cell1) + m.meshDensity(cell2)) / 2.0, 0.25);",
     "m.meshScalingDel2(iEdge) = 1.0 / libm::pow((m.meshDensity(cell1) + m.meshDensity(cell2)) / 4.0, 0.25);",
     D("meshScalingDel2", "SC JW")),
    (A, "m.meshScalingDel4(iEdge) = 1.0 / libm::pow((m.meshDensity(cell1) + m.meshDensity(cell2)) / 2.0, 0.75);",
     "m.meshScalingDel4(iEdge) = 1.0 / libm::pow((m.meshDensity(cell1) + m.meshDensity(cell2)) / 4.0, 0.75);",
     D("meshScalingDel4", "SC JW")),
    (A, "m.meshScalingRegionalEdge(iEdge) = 1.0 / libm::pow((m.meshDensity(cell1) + m.meshDensity(cell2)) / 2.0, 0.25);",
     "m.meshScalingRegionalEdge(iEdge) = 1.0 / libm::pow((m.meshDensity(cell1) + m.meshDensity(cell2)) / 4.0, 0.25);",
     D("meshScalingRegionalEdge", "SC JW")),
    (A, "      m.meshScalingRegionalCell(iCell) = 1.0 / libm::mvec_pow(m.meshDensity(iCell), 0.25);",
     "      m.meshScalingRegionalCell(iCell) = 1.0 / libm::mvec_pow(0.5 * m.meshDensity(iCell), 0.25);",
     D("meshScalingRegionalCell", "SC JW")),
    (A, "      if (m.cellsOnEdge(1, m.edgesOnCell(i, iCell)) != iCell) {", "      if (false) {",
     D("deformation_coef_c2 deformation_coef_s2 deformation_coef_cs", "MW SC JW")),
    (A, "      thetat[i] = libm::atan2(dy, dx) - pii / 2.;", "      thetat[i] = libm::atan2(dy, dx) + pii / 2.;",
     D("deformation_coef_c deformation_coef_s", "MW SC JW")),
    (A, "  fill_all(m.meshScalingDel2, 1.0);", "  fill_all(m.meshScalingDel2, 0.5);", D("meshScalingDel2", "MW SC")),
    (A, "  fill_all(m.meshScalingDel4, 1.0);", "  fill_all(m.meshScalingDel4, 0.5);", D("meshScalingDel4", "MW SC")),
    (A, "  fill_all(m.meshScalingRegionalCell, 1.0);", "  fill_all(m.meshScalingRegionalCell, 0.5);",
     D("meshScalingRegionalCell", "MW SC")),
    (A, "  fill_all(m.meshScalingRegionalEdge, 1.0);", "  fill_all(m.meshScalingRegionalEdge, 0.5);",
     D("meshScalingRegionalEdge", "MW SC")),
    (A, "      epssm_coeff_u[k] = min_coeff;", "      epssm_coeff_u[k] = max_coeff;", D("etp etm", "MW SC JW")),
    (A, "      epssm_coeff_w[k] = min_coeff;", "      epssm_coeff_w[k] = max_coeff;", D("ewp ewm", "MW SC JW")),
    (A, "    if (m.bdyMaskCell.data()[k] > nRelaxZone)", "    if (m.bdyMaskCell.data()[k] >= 0)", D("specZoneMaskCell", "MW SC JW")),
    (A, "    if (m.bdyMaskEdge.data()[k] > nRelaxZone)", "    if (m.bdyMaskEdge.data()[k] >= 0)", D("specZoneMaskEdge", "MW SC JW")),
    (A, "    if (m.bdyMaskVertex.data()[k] > nRelaxZone)", "    if (m.bdyMaskVertex.data()[k] >= 0)", D("specZoneMaskVertex", "MW SC JW")),
    (A, "  fill_all(m.nearestRelaxationCell, nCells + 1);", "  fill_all(m.nearestRelaxationCell, nCells);",
     D("nearestRelaxationCell", "MW SC JW")),
    (A, "      theta_m(k, iCell) = s.theta(k, iCell) * (1. + rvord * s.scalars(index_qv, k, iCell));",
     "      theta_m(k, iCell) = s.theta_base(k, iCell) * (1. + rvord * s.scalars(index_qv, k, iCell));", D("theta_m rtheta_p", "MW SC JW")),
    (A, "      rho_zz(k, iCell) = s.rho(k, iCell) / zz(k, iCell);", "      rho_zz(k, iCell) = s.rho(k, iCell) * (1.0 / zz(k, iCell));",
     D("rho_zz ru rho_p", "MW SC JW")),
    (A, "      d.ru(k, iEdge) = 0.5 * s.u(k, iEdge) * (rho_zz(k, cell1) + rho_zz(k, cell2));",
     "      d.ru(k, iEdge) = 0.5 * s.u(k, iEdge) * (rho_zz(k, cell1) + rho_zz(k, cell1));", D("ru", "MW SC JW")),
    (A, "(m.zb_cell(k, i, iCell) + std::copysign(1.0, flux) * m.zb3_cell(k, i, iCell))",
     "(m.zb_cell(k, i, iCell) + m.zb3_cell(k, i, iCell))", D("rw", "SC JW")),
    (A, "d.rho_p(k, iCell) = rho_zz(k, iCell) - s.rho_base(k, iCell);", "d.rho_p(k, iCell) = s.rho(k, iCell) - s.rho_base(k, iCell);",
     D("rho_p", "MW SC JW")),
    (A, "      d.exner(k, iCell) = vector ? libm::mvec_pow(xe, rcv) : libm::pow(xe, rcv);",
     "      d.exner(k, iCell) = libm::pow(xe, rcv);", D("exner pressure_p", "MW SC JW")),
    (A, "      d.exner_base(k, iCell) = vector ? libm::mvec_pow(xb, rcv) : libm::pow(xb, rcv);",
     "      d.exner_base(k, iCell) = libm::pow(xb, rcv);", D("exner_base pressure_base", "MW SC JW")),
    (A, "    for (int k = 1; k <= nVertLevels; ++k) d.rtheta_base(k, iCell) = s.theta_base(k, iCell) * s.rho_base(k, iCell);",
     "    for (int k = 1; k <= nVertLevels; ++k) d.rtheta_base(k, iCell) = s.rho_base(k, iCell) * s.theta_base(k, iCell) * 1.0000000000000002;",
     D("rtheta_base", "MW SC JW")),
    # SC: zb = zb3 = w = 0, zz = 1 and one areaTriangle, so divergence = rw = zb_cell = zb3_cell = 0 there; these
    # bugs (wrong array, missing sign) change them on every case
    (S, "      const real s = m.edgesOnCell_sign(i, iCell) * m.dvEdge(iEdge);", "      const real s = m.dvEdge(iEdge);",
     D("divergence", "MW SC JW")),
    (A, "m.invAreaTriangle(iVertex) = 1.0 / m.areaTriangle(iVertex);",
     "m.invAreaTriangle(iVertex) = 1.0 / m.kiteAreasOnVertex(1, iVertex);", D("invAreaTriangle", "MW SC JW")),
    (A, "      rho_zz(k, iCell) = s.rho(k, iCell) / zz(k, iCell);", "      rho_zz(k, iCell) = s.rho_base(k, iCell) / zz(k, iCell);",
     D("rho_zz rho_p", "MW SC JW")),
    (A, "    d.rw(1, iCell) = 0.0;", "    d.rw(1, iCell) = 1.0;", D("rw", "MW SC JW")),
    (A, "          m.zb_cell(k, i, iCell) = m.zb(k, side, e);", "          m.zb_cell(k, i, iCell) = m.zgrid(k, iCell);",
     D("zb_cell", "MW SC JW")),
    (A, "          m.zb3_cell(k, i, iCell) = m.zb3(k, side, e);", "          m.zb3_cell(k, i, iCell) = m.zgrid(k, iCell);",
     D("zb3_cell", "MW SC JW")),
    (W, "    std::vector<char> v(ref.c->begin(), ref.c->end());", "    std::vector<char> v(ref.c->rbegin(), ref.c->rend());",
     D("xtime", "MW SC JW")),
    ("tests/CMakeLists.txt", "  xtime\n  invAreaCell", "  invAreaCell", ["initdump_vars_MW", "initdump_vars_SC", "initdump_vars_JW"]),
    # --- review fixes of the review fixes: namelist values, integers, groups; atm_core_init; guards; writer; libm NaNs
    (NL, "    return bx == by;", "    return x == y;", ["unit_namelist"]),
    (NL, "    while (i < s.size() && (s[i] == ' ' || s[i] == '\\t' || s[i] == '\\r' || s[i] == '\\n')) {",
     "    while (i < s.size() && (s[i] == ' ' || s[i] == '\\t')) {", ["unit_namelist"]),
    (NL, "    if (errno == ERANGE || x < INT_MIN || x > INT_MAX) throw", "    if (false) throw", ["unit_namelist"]),
    (NL, "    if (v.text.size() > 512) throw", "    if (false) throw", ["unit_namelist"]),
    (NL, "      if (!groups.insert(record).second)", "      if (false)", ["unit_namelist"]),
    (NL, "    const auto at = [&] { return path + \":\" + std::to_string(name_line); };",
     "    const auto at = [&] { return path + \":\" + std::to_string(line); };", ["unit_namelist"]),
    (NL, "  if (c.config_microp_scheme == \"mp_kessler\") return {\"qv\", \"qc\", \"qr\"};",
     "  if (c.config_microp_scheme == \"mp_kessler\") return {\"qv\"};", ["unit_namelist", "err_scalars_missing"]),
    (NL, "  t(\"config_microp_scheme\", c.config_microp_scheme);", "", ["err_scalars_missing"]),
    (X, "    if (a.scalars != want) {", "    if (false) {", ["err_scalars_extra", "err_scalars_missing"]),
    (A, "  if (mo < 1 || mo > 12 || v(3) < 1 || v(3) > days[mo] || v(4) > 23 || v(5) > 59 || v(6) > 59)", "  if (false)",
     ["unit_init_guards"]),
    (A, "  const std::string stamp = atm_start_time_stamp(cfg);  // :238-239\n  atm_core_init_len_disp(m, cfg, log);                // :155-201",
     "  atm_core_init_len_disp(m, cfg, log);                // :155-201\n  const std::string stamp = atm_start_time_stamp(cfg);  // :238-239",
     ["unit_init_guards"]),
    (A, "  atm_core_init_len_disp(m, cfg, log);                // :155-201",
     "  { InitConfig c = cfg; atm_core_init_len_disp(m, c, log); }", ["unit_init_guards"]),
    (A, "  if (!m.bdyMaskCell.allocated()) throw", "  if (false) throw", ["unit_init_guards"]),
    (A, "  if (cfg.config_apply_lbcs)\n    throw Error(\"mpas_atm_run_compatibility: config_apply_lbcs = true (limited-area",
     "  if (false)\n    throw Error(\"mpas_atm_run_compatibility: config_apply_lbcs = true (limited-area", ["unit_init_guards"]),
    (A, "    if (m.meshDensity(iCell) != 1.0) return true;", "    if (false) return true;", ["unit_init_guards"]),
    (A, "    if (scaled_mesh && std::find(scaling.begin(), scaling.end(), n) != scaling.end()) continue;",
     "    if (false && std::find(scaling.begin(), scaling.end(), n) != scaling.end()) continue;", ["unit_init_guards"]),
    (A, "    if (n > 25) throw", "    if (n >= 25) throw", ["unit_init_guards"]),
    (V, "  if (N < 1) throw Error(\"mpas_legs: N < 1\");", "", ["unit_init_guards"]),
    (W, "  if (opt.nprocs < 1 || opt.nprocs > 9) throw", "  if (false) throw", ["unit_io_core"]),
    (W, "    std::remove(path.c_str());", "    (void)path;", ["unit_io_core"]),
    ("src/core/decomp.cpp",
     "  return static_cast<std::size_t>(j - 1) + static_cast<std::size_t>(dim1) * static_cast<std::size_t>(i - 1);",
     "  return static_cast<std::size_t>((j - 1) + dim1 * (i - 1));", ["unit_decomp"]),
    ("src/core/decomp.cpp", "    if (nEdgesOnCell[i - 1] < 0 || nEdgesOnCell[i - 1] > maxEdges)", "    if (false)", ["unit_decomp"]),
    ("src/core/decomp.cpp", "  if (nCellsGlobal < 0) throw", "  if (false) throw", ["unit_decomp"]),
    ("src/core/decomp.cpp", "  if (nCells < 0 || maxCells < 0 || nEdges < 0) throw", "  if (false) throw", ["unit_decomp"]),
    ("src/core/decomp.cpp", "  if (n2 < 0 || static_cast<std::size_t>(n2) > sorted.size())", "  if (false)", ["unit_decomp"]),
    ("src/core/decomp.cpp", "      innerDim > field.extent(0) || outerDim > field.extent(1))", "      false)", ["unit_decomp"]),
    ("src/core/decomp.cpp", "  if (innerDim < 0 || outerDim < 0 || indexSpaceDim < 0 || indexSpaceDim > indexToID.extent(0) ||",
     "  if (innerDim < 0 || outerDim < 0 || indexSpaceDim < 0 ||", ["unit_decomp"]),
    ("src/core/field.hpp", "      if (e != 0 && n > std::numeric_limits<std::size_t>::max() / static_cast<std::size_t>(e))",
     "      if (false)", ["unit_decomp"]),
    ("src/core/field.hpp", "    if (d < 0) throw Error(\"FArray: negative dimension index\");", "", ["unit_decomp"]),
    ("tests/compare_arrays.py", "            out.append((k, arr.dtype.str, arr.tobytes()))\n    return out",
     "            out.append((k, arr.dtype.str, arr.tobytes()))\n    return sorted(out)", ["compare_arrays_unit"]),
    # (an empty name would also fail as a missing variable: the mutant restores the old filter)
    ("tests/compare_arrays.py", "    names = a.vars.split(\",\")\n", "    names = [v for v in a.vars.split(\",\") if v]\n",
     ["compare_arrays_unit"]),
    ("tools/gen_namelist_table.py", "        if key in out:", "        if False:", ["namelist_gen_errors"]),
    ("tools/gen_namelist_table.py", "            if typ != \"character\" and quoted:", "            if False:", ["namelist_gen_errors"]),
    ("tools/gen_namelist_table.py", "        if not -2**31 <= int(text) <= 2**31 - 1:", "        if False:", ["namelist_gen_errors"]),
    ("tools/gen_namelist_table.py", "        if not math.isfinite(x):", "        if False:", ["namelist_gen_errors"]),
    ("tools/gen_namelist_table.py",
     "        t = re.match(r\"(?:[^!'\\\"]|'[^']*'|\\\"[^\\\"]*\\\")*\", line).group(0).strip()",
     "        t = line.split(\"!\")[0].strip()", ["namelist_gen_errors"]),
    (W, "      return s.substr(0, s.find('\\0'));", "      return s;", ["global_atts"]),
    (R, "  while (!s.empty() && s.back() == ' ') s.pop_back();\n  return s;\n}",
     "  while (!s.empty() && (s.back() == ' ' || s.back() == '\\0')) s.pop_back();\n  return s;\n}", ["global_atts"]),
    # NaN arguments: without the barrier GCC folds the caller's -x into the ports' x / x, x + y
    (L, "  asm(\"\" : \"+x\"(x));", "", ["unit_libm"]),
    (L, "1.0 ? __builtin_nan(\"\") : ::cppwrf::glibc::acos(x);", "1.0 ? -__builtin_nan(\"\") : ::cppwrf::glibc::acos(x);",
     ["unit_libm"]),
    # tests no mutant killed so far
    (R, "    if (!f.has_dim(\"nVertLevels\")) throw", "    if (false) throw", ["err_missing_nvl"]),
    (R, "    if (f.has_dim(\"nVertLevelsP1\") && f.dim_len", "    if (false && f.dim_len", ["err_nvlp1"]),
    (R, "      if (vi.shape[k] != kStrLen)", "      if (false)", ["err_strlen"]),
    (R, "  if (f.has_dim(\"TWO\") && f.dim_len(\"TWO\") != 2) throw", "  if (false) throw", ["err_two"]),
    # err_unknown_field: two guards (registry_var, then field_ref) give the same message, so no one-line change of
    # behaviour can reach the CLI (each guard is killed in unit_io_core); the message the user sees can be changed
    ("src/io/registry.cpp", "  if (!r) throw Error(\"variable '\" + name + \"' is not defined in MPAS-Atmosphere's Registry\");",
     "  if (!r) throw Error(\"variable '\" + name + \"' is not in the Registry\");", ["err_unknown_field", "unit_io_core"]),
]


# SC is planar without terrain (zb = zb3 = 0), has zz = 1, a uniform hexagonal mesh (one dcEdge, one areaTriangle),
# config_h_ScaleWithMesh = false and a horizontally uniform initial state for these terms: the mutants below are exact
# on SC (shown by mutation_check runs), so their SC tests are not required to fail; the other cases kill them.
EQUIVALENT_ON_SC = {
    ('src/core_atmosphere/atm_solve_diagnostics.cpp', 'divergence(k, iCell) = divergence(k, iCell) / m.areaCell(iCell);'):
        ['initdump_SC_divergence'],
    ('src/core_atmosphere/atm_solve_diagnostics.cpp', '(ke_edge(k, e3) + ke_edge(k, e2) + ke_edge(k, e1)) * r;'):
        ['initdump_SC_ke'],
    ('src/core_atmosphere/atm_solve_diagnostics.cpp', '(1. - ke_fact) * m.kiteAreasOnVertex(j, iVertex) * (ke_vertex(k, iVertex) * r);'):
        ['initdump_SC_ke'],
    ('src/core_atmosphere/atm_solve_diagnostics.cpp', 'pv_cell(k, iCell) + m.kiteAreasOnVertex(j, iVertex) * (pv_vertex(k, iVertex) * r);'):
        ['initdump_SC_gradPVn', 'initdump_SC_pv_cell'],
    ('src/core_atmosphere/atm_solve_diagnostics.cpp', '      const real s = m.edgesOnCell_sign(i, iCell) * m.dcEdge(iEdge);'):
        ['initdump_SC_divergence'],
    ('src/core_atmosphere/atm_solve_diagnostics.cpp', 'pv_cell(k, iCell) = 0.0;\n    const real r = m.invAreaCell(iCell);\n    for (int i = 1; i <= m.nEdgesOnCell(iCell); ++i) {\n      const int iVertex = m.verticesOnCell(1, iCell);'):
        ['initdump_SC_gradPVn'],
    ('src/core_atmosphere/atm_init.cpp', 'm.invDcEdge(iEdge) = 1.0 / std::sqrt(m.dcEdge(iEdge)) / std::sqrt(m.dcEdge(iEdge));'):
        ['initdump_SC_invDcEdge'],
    ('src/core_atmosphere/atm_init.cpp', 'm.invAreaTriangle(iVertex) = 1.0 / std::sqrt(m.areaTriangle(iVertex)) / std::sqrt(m.areaTriangle(iVertex));'):
        ['initdump_SC_invAreaTriangle'],
    ('src/core_atmosphere/atm_init.cpp', '          m.zb3_cell(k, i, iCell) = m.zb(k, side, e);'):
        ['initdump_SC_zb3_cell'],
    ('src/core_atmosphere/atm_init.cpp', '        d.rw(k, iCell) = d.rw(k, iCell) - '):
        ['initdump_SC_rw'],
    ('src/core_atmosphere/atm_init.cpp', '          m.zb_cell(k, i, iCell) = m.zb(k, 1, e);'):
        ['initdump_SC_zb_cell'],
    ('src/core_atmosphere/atm_init.cpp', '          m.zb3_cell(k, i, iCell) = m.zb3(k, 1, e);'):
        ['initdump_SC_zb3_cell'],
    ('src/core_atmosphere/atm_init.cpp', 'm.meshScalingDel2(iEdge) = 1.0 / libm::pow((m.meshDensity(cell1) + m.meshDensity(cell2)) / 4.0, 0.25);'):
        ['initdump_SC_meshScalingDel2'],
    ('src/core_atmosphere/atm_init.cpp', 'm.meshScalingDel4(iEdge) = 1.0 / libm::pow((m.meshDensity(cell1) + m.meshDensity(cell2)) / 4.0, 0.75);'):
        ['initdump_SC_meshScalingDel4'],
    ('src/core_atmosphere/atm_init.cpp', 'm.meshScalingRegionalEdge(iEdge) = 1.0 / libm::pow((m.meshDensity(cell1) + m.meshDensity(cell2)) / 4.0, 0.25);'):
        ['initdump_SC_meshScalingRegionalEdge'],
    ('src/core_atmosphere/atm_init.cpp', '      m.meshScalingRegionalCell(iCell) = 1.0 / libm::mvec_pow(0.5 * m.meshDensity(iCell), 0.25);'):
        ['initdump_SC_meshScalingRegionalCell'],
    ('src/core_atmosphere/atm_init.cpp', '      rho_zz(k, iCell) = s.rho(k, iCell) * (1.0 / zz(k, iCell));'):
        ['initdump_SC_rho_p', 'initdump_SC_rho_zz', 'initdump_SC_ru'],
    ('src/core_atmosphere/atm_init.cpp', '(m.zb_cell(k, i, iCell) + m.zb3_cell(k, i, iCell))'):
        ['initdump_SC_rw'],
    ('src/core_atmosphere/atm_init.cpp', 'd.rho_p(k, iCell) = s.rho(k, iCell) - s.rho_base(k, iCell);'):
        ['initdump_SC_rho_p'],
}
MUTATIONS = [(p, o, n, [t for t in ts if t not in EQUIVALENT_ON_SC.get((p, n), [])]) for p, o, n, ts in MUTATIONS]


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
