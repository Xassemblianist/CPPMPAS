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
RS = "src/io/mpas_restart.cpp"


def RST0(fields):
    """step-0 tests of the restart run (tests/CMakeLists.txt): the state written is the restart stream's"""
    return [f"dycore_restart0_MW_{f}" for f in fields.split()]

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
    (W, "          f.put_att_text(n, k, v);\n        }",
     "          f.put_att_text(n, k == \"long_name\" ? \"description\" : k, v);\n        }",
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
    # the reader's cut of parent_id at StrKIND (domain % parent_id) is not visible in the output since the writer
    # cuts every text attribute at StrKIND as MPAS_io_put_att_text does: that cut is the one mutated
    (W, "      s = s.substr(0, 512);\n", "", ["global_atts"]),
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
    (NL, "  t(\"config_start_time\", c.config_start_time);", "",
     ["init_step0_MW", "init_step0_SC", "init_step0_JW"] + [f"dycore_MW_step{n}_xtime" for n in range(11)] +
     ["dycore_MW_restart_xtime", "dycore_MW_output30_xtime"]),
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
     ["readrestart_MW", "readrestart_SC", "readrestart_JW"] + [f"dycore_restart0_MW_{f}" for f in "u w rho theta".split()]),
    ("tests/compare_arrays.py", "            out.append((k, arr.dtype.str, arr.tobytes()))", "            out.append((k, \"x\", str(a)))",
     ["compare_arrays_unit"]),
    ("tests/compare_arrays.py", "            if attributes(x) != attributes(y):", "            if False:", ["compare_arrays_unit"]),
    (L, "double cos(double x) { return ::cppwrf::glibc::cos(arg(x)); }",
     "double cos(double x) { return ::cppwrf::glibc::mvec::cos(arg(x)); }", ["unit_libm"]),
    (L, "    *sinx = g::copysign(g::do_cos(a, da), x);", "    *sinx = g::copysign(g::do_cos(y, hp1), x);", ["unit_libm"]),
    (L, "1.0 ? __builtin_nan(\"\") : ::cppwrf::glibc::asin(x);", "1.0 ? -__builtin_nan(\"\") : ::cppwrf::glibc::asin(x);",
     ["unit_libm"]),
    # --- every initdump field (atm_mpas_init_block at t = 0) against mpas842-dbl-mw-003 / sc-002 / jw-002
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
     D("pv_cell gradPVn", "SC JW") + D("pv_cell", "MW")),  # MW gradPVn: neighbouring cells change alike
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
    (RC, "            -((uReconstructX(k, iCell) * clon + uReconstructY(k, iCell) * slon) * slat) + uReconstructZ(k, iCell) * clat;",
     "            -(uReconstructX(k, iCell) * clon * slat + uReconstructY(k, iCell) * slon * slat) + uReconstructZ(k, iCell) * clat;",
     D("uReconstructMeridional", "JW")),
    (RC, "      libm::sincos(m.lonCell(iCell), &slon, &clon);", "      libm::sincos(m.lonCell(iCell), &clon, &slon);",
     D("uReconstructZonal uReconstructMeridional", "JW")),
    # the init-block arrays the restart and step files do not hold (initdump only)
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
    # --- the two routines on arrays that hold values (unit_solve_diag): each zeroing of an accumulation
    (S, "    for (int k = 1; k <= nVertLevels; ++k) vorticity(k, iVertex) = 0.0;\n", "", ["unit_solve_diag"]),
    (S, "    for (int k = 1; k <= nVertLevels; ++k) divergence(k, iCell) = 0.0;\n", "", ["unit_solve_diag"]),
    (S, "    for (int k = 1; k <= nVertLevels; ++k) ke(k, iCell) = 0.0;\n", "", ["unit_solve_diag"]),
    (S, "      for (int k = 1; k <= nVertLevels; ++k) v(k, iEdge) = 0.0;\n", "", ["unit_solve_diag"]),
    (S, "    for (int k = 1; k <= nVertLevels; ++k) pv_cell(k, iCell) = 0.0;\n", "", ["unit_solve_diag"]),
    (RC, "      uReconstructX(k, iCell) = 0.0;\n", "", ["unit_solve_diag"]),
    (RC, "      uReconstructY(k, iCell) = 0.0;\n", "", ["unit_solve_diag"]),
    (RC, "      uReconstructZ(k, iCell) = 0.0;\n", "", ["unit_solve_diag"]),
    # rk_step (TI:7307-7310): always recomputing v, never, inverted, the absent case, the range check
    (S, "    if (*rk_step != 3) reconstruct_v = false;", "    if (*rk_step == 0) reconstruct_v = false;", ["unit_solve_diag"]),
    (S, "    if (*rk_step != 3) reconstruct_v = false;", "    if (*rk_step != 1) reconstruct_v = false;", ["unit_solve_diag"]),
    (S, "  if (rk_step.has_value()) {", "  if (true) {", ["unit_solve_diag"]),
    (S, "  if (rk_step && (*rk_step < 1 || *rk_step > 3))", "  if (rk_step && (*rk_step < 0 || *rk_step > 3))",
     ["unit_solve_diag"]),
    (S, "  if (rk_step && (*rk_step < 1 || *rk_step > 3))", "  if (rk_step && (*rk_step < 1 || *rk_step > 4))",
     ["unit_solve_diag"]),
    # one thread's ranges, extents, kiteForCell, vertexDegree, APVM, nVertLevels
    (S, "  if (cellStart != 1 || cellEnd != m.nCells ||", "  if (cellEnd != m.nCells ||", ["unit_solve_diag"]),
    (S, "vertexEnd != m.nVertices || edgeStart != 1 ||", "vertexEnd != m.nVertices ||", ["unit_solve_diag"]),
    (S, "  if (a.rank() != 2 || a.extent(0) != nL || a.extent(1) != n1)\n    fail(", "  if (a.rank() != 2 || a.extent(0) != nL)\n    fail(",
     ["unit_solve_diag"]),
    (S, "      if (m.kiteForCell(i, iCell) < 1 || m.kiteForCell(i, iCell) > vertexDegree)",
     "      if (m.kiteForCell(i, iCell) < 1)", ["unit_solve_diag"]),
    (S, "      if (m.kiteForCell(i, iCell) < 1 || m.kiteForCell(i, iCell) > vertexDegree)",
     "      if (m.kiteForCell(i, iCell) > vertexDegree)", ["unit_solve_diag"]),
    (S, "  if (vertexDegree != 3) fail(", "  if (false) fail(", ["unit_solve_diag"]),
    (S, "  if (!(config_apvm_upwinding > 0.0))", "  if (config_apvm_upwinding == 0.0)", ["unit_solve_diag"]),
    (S, "  if (nL < 1) throw Error(\"allocate_solve_diag_work", "  if (false) throw Error(\"allocate_solve_diag_work",
     ["unit_solve_diag"]),
    (S, "  if (nL < 1) throw Error(\"allocate_solve_diag_fields", "  if (false) throw Error(\"allocate_solve_diag_fields",
     ["unit_solve_diag"]),
    (S, "  w.ke_vertex.allocate({nL, m.nVertices + 1}, 0.0);", "  w.ke_vertex.allocate({nL, m.nVertices + 2}, 0.0);",
     ["unit_solve_diag"]),
    (S, "  for (auto* a : {&d.vorticity, &d.pv_vertex}) a->allocate({nL, nV}, 0.0);",
     "  for (auto* a : {&d.vorticity, &d.pv_vertex}) a->allocate({nL, nV + 1}, 0.0);", ["unit_solve_diag"]),
    # mpas_reconstruct_2d: Fortran's -a*b is -(a*b) (no negation of a NaN); sincos, not sin / cos; its checks
    (RC, "-(uReconstructX(k, iCell) * slon) + uReconstructY(k, iCell) * clon;",
     "-uReconstructX(k, iCell) * slon + uReconstructY(k, iCell) * clon;", ["unit_solve_diag"]),
    (RC, "-((uReconstructX(k, iCell) * clon + uReconstructY(k, iCell) * slon) * slat) + uReconstructZ(k, iCell) * clat;",
     "-(uReconstructX(k, iCell) * clon + uReconstructY(k, iCell) * slon) * slat + uReconstructZ(k, iCell) * clat;",
     ["unit_solve_diag"]),
    # (sincos -> sin / cos of latCell is not here: on |lat| < pi/2 no value is known where they differ, see
    # tests/test_solve_diag.cpp)
    (RC, "      libm::sincos(m.lonCell(iCell), &slon, &clon);",
     "      slon = libm::sin(m.lonCell(iCell)), clon = libm::cos(m.lonCell(iCell));", ["unit_solve_diag"]),
    (RC, "  if (a.rank() != 2 || a.extent(0) != nL || a.extent(1) != n1)\n    throw",
     "  if (a.rank() != 2 || a.extent(0) != nL)\n    throw", ["unit_solve_diag"]),
    (RC, "coeffs_reconstruct.extent(2) < nCells)\n    throw Error(",
     "coeffs_reconstruct.extent(2) < nCells)\n    (void)0;\n  if (false) throw Error(", ["unit_solve_diag"]),
    # the harness itself: a test that was not run is not a kill
    ("tests/mutation_check.py", "    return status is not None and status != \"Passed\" and \"Not Run\" not in status",
     "    return status is not None and status != \"Passed\"", ["mutation_check_unit"]),
    # --unverified on the CLI (sc_init.nc: nVertLevels = 3) and the exact verified lists
    (X, "    else if (s == \"--unverified\") a.unverified = true;", "    else if (s == \"--unverified\") a.unverified = false;",
     ["init_unverified_SC"]),
    (X, "      if (is_computed && !a.unverified && std::find(verified.begin(), verified.end(), f) == verified.end())",
     "      if (false)", ["err_unverified"]),
    (A, "  const std::vector<std::string> pow_tail = {\"exner\", \"exner_base\", \"pressure_p\", \"pressure_base\"};",
     "  const std::vector<std::string> pow_tail = {\"exner\", \"pressure_p\", \"pressure_base\"};", ["unit_init_guards"]),
]


# ---------------------------------------------------------------------------------------------------------------
# the dycore (src/core_atmosphere/atm_dycore.cpp, atm_dissipation.cpp) against the MW step files, its checks and the
# run command
DY, DS = "src/core_atmosphere/atm_dycore.cpp", "src/core_atmosphere/atm_dissipation.cpp"


def ST(steps, fields="u w rho_zz theta_m rho theta pressure_p"):
    """dycore_MW_step<n>_<field> tests (tests/CMakeLists.txt)"""
    return [f"dycore_MW_step{n}_{f}" for n in steps for f in fields.split()]


S1 = ST([1])
ALLDYN = ST(range(1, 11)) + ["dycore_MW_mpas_compare"]
LONG_R = "u w rho_zz theta_m pressure_p rho theta exner ru rw rtheta_p rho_p".split()
LONG_O = "u w rho theta divergence vorticity ke uReconstructZonal uReconstructMeridional".split()
LONG = [f"dycore_MW_restart_{f}" for f in LONG_R] + [f"dycore_MW_output30_{f}" for f in LONG_O]
MUTATIONS += [
    # atm_rk_integration_setup, moist coefficients, physics_get_tend, atm_compute_vert_imp_coefs
    (DY, "      f.u2(k, iEdge) = s.u(k, iEdge);", "      f.u2(k, iEdge) = 0.999 * s.u(k, iEdge);", S1),
    (DY, "      f.rw_save(k, iCell) = d.rw(k, iCell);", "      f.rw_save(k, iCell) = 0.0;", S1),
    (DY, "      f.cqu(k, iEdge) = 1.0 / (1.0 + qtotal);", "      f.cqu(k, iEdge) = 1.0 / (1.01 + qtotal);", S1),
    (DY, "      f.cqw(k, iCell) = 1.0 / (1.0 + qtotal);", "      f.cqw(k, iCell) = 1.0 / (1.01 + qtotal);", S1),
    (DY, "      f.tend_rtheta_physics(k, i) = f.tend_rtheta_physics(k, i) + tend_th;",
     "      f.tend_rtheta_physics(k, i) = f.tend_rtheta_physics(k, i) + tend_th + 1.0e-9;", S1),
    (DY, "      alpha_tri(k, iCell) = 1. / (1.0 + (dts * dts) *", "      alpha_tri(k, iCell) = 1. / (1.0 + (dts * 2.0) *", S1),
    (DY, "      gamma_tri(k, iCell) = (dts * dts) * c_tri[k] * alpha_tri(k, iCell);",
     "      gamma_tri(k, iCell) = (dts * dts) * c_tri[k] * alpha_tri(k - 1, iCell);", S1),
    # (etp and ewp are 0.5 * (1 + epssm) on every level of MW: a mutant of their level index is exact there)
    (DY, "etp(k - 1) * rdzw(k - 1) * zz(k - 1, iCell))", "etp(k - 1) * rdzw(k - 1) * zz(k, iCell))", S1),
    (DY, "cofwt(k, iCell) * coftz(k + 1, iCell) * rdzw(k);", "cofwt(k, iCell) * coftz(k, iCell) * rdzw(k);", S1),
    (DY, "cofwr(k, iCell) * cofrz(k - 1) - cofwt(k - 1, iCell) * coftz(k - 1, iCell) * rdzw(k - 1);",
     "cofwr(k, iCell) * cofrz(k - 1) - cofwt(k, iCell) * coftz(k - 1, iCell) * rdzw(k - 1);", S1),
    (DY, "* rb(k, iCell) / (1. + qtotal) * p(k, iCell) /", "* rb(k, iCell) / (1.1 + qtotal) * p(k, iCell) /", S1),
    (DY, "      coftz(k, iCell) = (fzm(k) * t(k, iCell) + fzp(k) * t(k - 1, iCell));",
     "      coftz(k, iCell) = (fzm(k) * t(k, iCell) + fzp(k) * t(k, iCell));", S1),
    (DY, "cofwr(k, iCell) = .5 * gravity *", "cofwr(k, iCell) = .25 * gravity *", S1),
    # atm_compute_dyn_tend
    (DY, "      for (int k = 1; k <= nVertLevels; ++k) h_divergence(k, iCell) = h_divergence(k, iCell) + edge_sign * ru(k, iEdge);",
     "      for (int k = 1; k <= nVertLevels; ++k) h_divergence(k, iCell) = h_divergence(k, iCell) + edge_sign * u(k, iEdge);", S1),
    (DY, "        tend_rho(k, iCell) = -h_divergence(k, iCell) - rdzw(k) * (rw(k + 1, iCell) - rw(k, iCell)) +",
     "        tend_rho(k, iCell) = -h_divergence(k, iCell) - rdzw(k) * (rw(k + 1, iCell) + rw(k, iCell)) +", S1),
    (DY, "rr_save(k, iCell) * (1. + f.qtot(k, iCell))", "rr_save(k, iCell) * (1.5 + f.qtot(k, iCell))", S1),
    (DY, "                                 0.5 * m.zxu(k, iEdge) * (dpdz(k, cell1) + dpdz(k, cell2))));",
     "                                 0.25 * m.zxu(k, iEdge) * (dpdz(k, cell1) + dpdz(k, cell2))));", S1),
    (DY, "    int k = 2;\n    wduz[k] = 0.5 *", "    int k = 2;\n    wduz[k] = 0.25 *", S1),
    (DY, "    for (k = 1; k <= nVertLevels; ++k) tend_u(k, iEdge) = -(rdzw(k) * (wduz[k + 1] - wduz[k]));",
     "    for (k = 1; k <= nVertLevels; ++k) tend_u(k, iEdge) = -(rdzw(k) * (wduz[k + 1] + wduz[k]));", S1),
    (DY, "        const real workpv = 0.5 * (pv_edge(k, iEdge) + pv_edge(k, eoe));",
     "        const real workpv = 0.5 * (pv_edge(k, iEdge) + pv_edge(k, iEdge));", S1),
    (DY, "(q[k] - (ke(k, cell2) - ke(k, cell1)) * m.invDcEdge(iEdge))", "(q[k] - (ke(k, cell2) + ke(k, cell1)) * m.invDcEdge(iEdge))", S1),
    (DY, "      tend_u(k, iEdge) = tend_u(k, iEdge) + tend_u_euler(k, iEdge) + f.tend_ru_physics(k, iEdge);",
     "      tend_u(k, iEdge) = tend_u(k, iEdge) + f.tend_ru_physics(k, iEdge);", S1),
    (DY, "          const real scalar_weight = m.adv_coefs(j, iEdge) + sign1(ru_edge_w[k]) * m.adv_coefs_3rd(j, iEdge);",
     "          const real scalar_weight = m.adv_coefs(j, iEdge) - sign1(ru_edge_w[k]) * m.adv_coefs_3rd(j, iEdge);", S1),
    (DY, "      tend_w(k, iCell) = tend_w(k, iCell) * m.invAreaCell(iCell) - rdzu(k) * (wdwz[k + 1] - wdwz[k]);",
     "      tend_w(k, iCell) = tend_w(k, iCell) * m.invAreaCell(iCell) + rdzu(k) * (wdwz[k + 1] - wdwz[k]);", S1),
    (DY, "            tend_w_euler(k, iCell) - f.cqw(k, iCell) * (rdzu(k) * (pp(k, iCell) - pp(k - 1, iCell)) -",
     "            tend_w_euler(k, iCell) - f.cqw(k, iCell) * (rdzu(k) * (pp(k, iCell) - pp(k, iCell)) -", S1),
    (DY, "    for (int k = 2; k <= nVertLevels; ++k) tend_w(k, iCell) = tend_w(k, iCell) + tend_w_euler(k, iCell);",
     "    for (int k = 2; k <= nVertLevels; ++k) tend_w(k, iCell) = tend_w(k, iCell);", S1),
    (DY, "          const real scalar_weight = m.adv_coefs(j, iEdge) + sign1(ru(k, iEdge)) * m.adv_coefs_3rd(j, iEdge);",
     "          const real scalar_weight = m.adv_coefs(j, iEdge) - sign1(ru(k, iEdge)) * m.adv_coefs_3rd(j, iEdge);", S1),
    (DY, "(ru_save(k, iEdge) - ru(k, iEdge)) * 0.5 *", "(ru_save(k, iEdge) - ru(k, iEdge)) * 0.25 *", S1),
    (DY, "                      coef_3rd_order);", "                      2.0 * coef_3rd_order);", S1),
    (DY, "    wdtz[k] = rw_save(k, iCell) * (fzm(k) * theta_m(k, iCell) + fzp(k) * theta_m(k - 1, iCell));",
     "    wdtz[k] = rw(k, iCell) * (fzm(k) * theta_m(k, iCell) + fzp(k) * theta_m(k - 1, iCell));", S1),
    (DY, "      tend_theta(k, iCell) = tend_theta(k, iCell) + tend_theta_euler(k, iCell) + f.tend_rtheta_physics(k, iCell);",
     "      tend_theta(k, iCell) = tend_theta(k, iCell) + f.tend_rtheta_physics(k, iCell);", S1),
    # dissipation
    (DS, "      const real kdiffu = 0.5 * (eddy_visc_horz(k, cell1) + eddy_visc_horz(k, cell2));",
     "      const real kdiffu = 0.25 * (eddy_visc_horz(k, cell1) + eddy_visc_horz(k, cell2));", S1),
    (DS, "                                      (u(k, iEdge) - u(k - 1, iEdge)) / (z0 - zm)) /",
     "                                      (u(k, iEdge) - u(k - 1, iEdge)) / (z0 - zp)) /", ST([2])),
    (DS, "        w_turb_flux = w_turb_flux * m.meshScalingDel2(iEdge) * 0.25 *",
     "        w_turb_flux = w_turb_flux * m.meshScalingDel2(iEdge) * 0.5 *", S1),
    (DS, "                                                               (w(k, iCell) - w(k - 1, iCell)) * m.rdzw(k - 1)) *",
     "                                                               (w(k, iCell) - w(k - 1, iCell)) * m.rdzw(k)) *", S1),
    (DS, "            theta_turb_flux * 0.5 * (eddy_visc_horz(k, cell1) + eddy_visc_horz(k, cell2)) * pr_scale;  // DM:1295",
     "            theta_turb_flux * 0.25 * (eddy_visc_horz(k, cell1) + eddy_visc_horz(k, cell2)) * pr_scale;  // DM:1295", S1),
    (DS, "                                          (theta_m(k, iCell) - theta_m(k - 1, iCell)) / (z0 - zm)) /",
     "                                          (theta_m(k, iCell) - theta_m(k - 1, iCell)) / (z0 - zp)) /", S1),
    # atm_set_smlstep_pert_variables: the sign is that of u_tend(k), not of the flux (TI:3629)
    (DY, "                           (m.zb_cell(k, i, iCell) + sign1(u_tend(k, iEdge)) * m.zb3_cell(k, i, iCell)) * flux;",
     "                           (m.zb_cell(k, i, iCell) + sign1(flux) * m.zb3_cell(k, i, iCell)) * flux;", S1),
    (DY, "      w_tend(k, iCell) = (m.fzm(k) * m.zz(k, iCell) + m.fzp(k) * m.zz(k - 1, iCell)) * w_tend(k, iCell);",
     "      w_tend(k, iCell) = (m.fzm(k) * m.zz(k, iCell) + m.fzp(k) * m.zz(k, iCell)) * w_tend(k, iCell);", S1),
    # atm_advance_acoustic_step, atm_divergence_damping_3d
    (DY, "        pgrad = cqu(k, iEdge) * 0.5 * c2 *", "        pgrad = cqu(k, iEdge) * 0.25 * c2 *", S1),
    (DY, "        pgrad = pgrad + 0.5 * m.zxu(k, iEdge) * gravity", "        pgrad = pgrad + 0.25 * m.zxu(k, iEdge) * gravity", S1),
    (DY, "      for (int k = 1; k <= nVertLevels; ++k) rtheta_pp_old(k, iCell) = rtheta_pp(k, iCell);",
     "      for (int k = 1; k <= nVertLevels; ++k) rtheta_pp_old(k, iCell) = 0.0;", S1),
    (DY, "      rs[k] = rho_pp(k, iCell) + dts * tend_rho(k, iCell) + rs[k] -", "      rs[k] = rho_pp(k, iCell) + rs[k] -",
     ALLDYN + LONG),
    (DY, "        ts[k] = ts[k] - flux * 0.5 * (theta_m(k, cell2) + theta_m(k, cell1));",
     "        ts[k] = ts[k] - flux * 0.5 * (theta_m(k, cell2) + theta_m(k, cell2));", S1),
    (DY, "                 cofwt(k - 1, iCell) * (etp(k - 1) * ts[k - 1]", "                 cofwt(k, iCell) * (etp(k - 1) * ts[k - 1]", S1),
    (DY, "      rw_p(k, iCell) = rw_p(k, iCell) - gamma_tri(k, iCell) * rw_p(k + 1, iCell);",
     "      rw_p(k, iCell) = rw_p(k, iCell) + gamma_tri(k, iCell) * rw_p(k + 1, iCell);", S1),
    (DY, "                           (1.0 + dts * dss(k, iCell)) -", "                           (1.0 + 2.0 * dts * dss(k, iCell)) -", S1),
    (DY, "      rho_pp(k, iCell) = rs[k] - dts * cofrz(k) * (ewp(k + 1) * rw_p(k + 1, iCell) - ewp(k) * rw_p(k, iCell));",
     "      rho_pp(k, iCell) = rs[k] - dts * cofrz(k) * (ewm(k + 1) * rw_p(k + 1, iCell) - ewm(k) * rw_p(k, iCell));", S1),
    (DY, "  const real coef_divdamp = 2.0 * cfg.config_smdiv * cfg.config_len_disp * rdts;",
     "  const real coef_divdamp = 4.0 * cfg.config_smdiv * cfg.config_len_disp * rdts;", S1),
    # atm_recover_large_step_variables
    (DY, "      rho_zz(k, iCell) = rho_p(k, iCell) + rho_base(k, iCell);", "      rho_zz(k, iCell) = rho_p(k, iCell) + 1.000001 * rho_base(k, iCell);", S1),
    (DY, "      w(k, iCell) = rw(k, iCell) / (fzm(k) * zz(k, iCell) + fzp(k) * zz(k - 1, iCell));",
     "      w(k, iCell) = rw(k, iCell) / (fzm(k) * zz(k, iCell) + fzp(k) * zz(k, iCell));", S1),
    (DY, "        exner(k, iCell) = vector ? libm::mvec_pow(x, rcv) : libm::pow(x, rcv);",
     "        exner(k, iCell) = vector ? libm::mvec_pow(x, 1.0001 * rcv) : libm::pow(x, rcv);",
     ["dycore_MW_step1_pressure_p"] + ST([2])),
    # the rk_step 3 loop's level pairs call libmvec's pow: scalar pow on every level differs from step 1's pressure_p
    (DY, "        exner(k, iCell) = vector ? libm::mvec_pow(x, rcv) : libm::pow(x, rcv);",
     "        exner(k, iCell) = libm::pow(x, rcv);", ["dycore_MW_step1_pressure_p"] + ST(range(2, 11))),
    (DY, "(exner(k, iCell) * rtheta_p(k, iCell) + rtheta_base(k, iCell) * (exner(k, iCell) - exner_base(k, iCell)));",
     "(exner(k, iCell) * rtheta_p(k, iCell) + rtheta_base(k, iCell) * (exner(k, iCell) - exner(k, iCell)));",
     ["dycore_MW_step1_pressure_p"] + ST([2])),
    (DY, "      u(k, iEdge) = 2. * ru(k, iEdge) / (rho_zz(k, cell1) + rho_zz(k, cell2));",
     "      u(k, iEdge) = 2. * ru(k, iEdge) / (rho_zz(k, cell1) + rho_zz(k, cell1));", S1),
    (DY, "      real flux = (cf1 * ru(1, iEdge) + cf2 * ru(2, iEdge) + cf3 * ru(3, iEdge));",
     "      real flux = (cf1 * ru(1, iEdge) + cf2 * ru(2, iEdge));", S1),
    (DY, "    w(1, iCell) = w(1, iCell) / (cf1 * rho_zz(1, iCell) + cf2 * rho_zz(2, iCell) + cf3 * rho_zz(3, iCell));",
     "    w(1, iCell) = w(1, iCell) / (cf1 * rho_zz(1, iCell) + cf2 * rho_zz(2, iCell));", S1),
    (DY, "      theta_m(k, iCell) = (rtheta_p(k, iCell) + rtheta_base(k, iCell)) / rho_zz(k, iCell);\n      }",
     "      theta_m(k, iCell) = (rtheta_p(k, iCell) + rtheta_base(k, iCell)) / rho_zz(k, iCell) * 1.000001;\n      }", S1),
    # atm_srk3: step constants, solve diagnostics per rk_step, the scalar clamp
    (DY, "  rk_sub_timestep[1] = dt_dynamics / static_cast<real>(number_of_sub_steps);",
     "  rk_sub_timestep[1] = dt_dynamics / static_cast<real>(number_of_sub_steps + 1);", S1),
    (DY, "  number_sub_steps[3] = number_of_sub_steps;", "  number_sub_steps[3] = number_of_sub_steps - 1;", S1),
    (DY, "    if (f.scalars2.data()[i] < 0.0) f.scalars2.data()[i] = 0.0;",
     "    if (f.scalars2.data()[i] <= 0.0) f.scalars2.data()[i] = 1.0e-300;",
     ST(range(1, 11), "qv") + ["dycore_MW_restart_qv", "dycore_MW_output30_qv", "dycore_restart_MW_output60_qv"]),
    # output diagnostics, Time, xtime, the step files
    (DY, "      s.rho(k, iCell) = s.rho_zz(k, iCell) * m.zz(k, iCell);", "      s.rho(k, iCell) = s.rho_zz(k, iCell);",
     ST(range(11), "rho")),
    (DY, "      s.theta(k, iCell) = s.theta_m(k, iCell) / (1. + rvord * s.scalars(index_qv, k, iCell));",
     "      s.theta(k, iCell) = s.theta_m(k, iCell) / (1.5 + rvord * s.scalars(index_qv, k, iCell));", ST(range(11), "theta")),
    (DY, "  f.Time2 = f.Time1 + dt;", "  f.Time2 = f.Time1 + dt + 1.0;", ST(range(1, 11), "Time") + ["dycore_MW_restart_Time"]),
    ("src/core_atmosphere/atm_dycore.hpp", "  real Time1 = 0.0, Time2 = 0.0;", "  real Time1 = 1.0, Time2 = 0.0;",
     ST(range(11), "Time")),
    (X, "      atm_timestep(dcfg, in.mesh, in.state, in.diag, f, dcfg.config_dt, add_seconds(start, n * dt_s));",
     "      atm_timestep(dcfg, in.mesh, in.state, in.diag, f, dcfg.config_dt, add_seconds(start, (n + 1) * dt_s));",
     ST([1], "xtime u w") + ["dycore_restart_MW_output60_xtime"]),
    (X, "    write_step(0);", "", ST([0], "xtime Time u w rho_zz theta_m rho theta pressure_p qv")),
    (X, "    for (const auto& sn : in.state.scalar_names) fields.push_back(sn);", "", ST(range(11), "qv")),
    ("src/io/mpas_writer.cpp", '          f.put_att_text(n, k, "seconds since " + t);', "          f.put_att_text(n, k, v);",
     ST(range(11), "Time")),
    # the configuration checks, the CLI
    (DY, '  if (c.config_horiz_mixing != "2d_fixed") no(', '  if (c.config_horiz_mixing != "2d_smagorinsky") no(',
     ["dycore_run_MW", "dycore_run_long_MW", "dycore_restart0_run_MW", "dycore_restart_run_long_MW", "err_dycore_JW", "unit_dycore"]),
    (DY, "  if (c.config_scalar_advection) no(", "  if (false) no(", ["err_dycore_SC", "unit_dycore"]),
    (DY, "    if (m.fEdge(iEdge) != 0.0) no(", "    if (false) no(", ["err_dycore_fEdge", "unit_dycore"]),
    (DY, "  if (m.nVertLevels < 4 || m.nVertLevels % 2 != 0) no(", "  if (false) no(", ["err_dycore_odd_levels", "unit_dycore"]),
    (DY, "  if (c.config_time_integration_order != 2)", "  if (c.config_time_integration_order > 3)", ["unit_dycore"]),
    (DY, "  if (!(c.config_v_mom_eddy_visc2 > 0.0)) no(", "  if (false) no(", ["unit_dycore"]),
    (DY, "  if (!off(c.config_microp_scheme)) no(", "  if (false) no(", ["unit_dycore"]),
    (DY, "  if (!(c.config_dt > 0.0) || c.config_dt != std::floor(c.config_dt) || c.config_dt > 86400.0)",
     "  if (!(c.config_dt > 0.0) || c.config_dt > 86400.0)", ["unit_dycore"]),
    (DY, "  c.config_len_disp = init_cfg.config_len_disp;", '  c.config_len_disp = r("config_len_disp");', ["unit_dycore"]),
    (DY, "  auto leap = [](int y) { return (y % 4 == 0 && y % 100 != 0) || y % 400 == 0; };",
     "  auto leap = [](int y) { return y % 4 == 0; };", ["unit_dycore"]),
    (DY, "    if (++D > mdays(Y, M)) {", "    if (++D >= mdays(Y, M)) {", ["unit_dycore"]),
    (X, "  if (used != v.size() || n < 0 || n > 1000000) throw Error(opt +", "  if (n < 0 || n > 1000000) throw Error(opt +",
     ["err_dycore_steps"]),
    (DY, "  mpas_reconstruct_2d(m, f.u2, d.uReconstructX,", "  if (false) mpas_reconstruct_2d(m, f.u2, d.uReconstructX,",
     ["dycore_MW_output30_uReconstructZonal", "dycore_MW_output30_uReconstructMeridional"]),
    (X, "      if (n > a.steps) throw Error(\"--write-at \"", "      if (false) throw Error(\"--write-at \"", ["err_dycore_write_at"]),
    (X, "    for (const auto& w : a.write_fields) (void)field_ref(in, w);", "", ["err_dycore_write_fields"]),
    (X, "    if (a.pos.size() != 2 || a.steps < 0 || a.namelist.empty() || !a.fields.empty())",
     "    if (a.pos.size() != 2 || a.namelist.empty() || !a.fields.empty())", ["err_dycore_no_steps"]),
    # the restart run (src/io/mpas_restart.cpp, the restart branch of atm_mpas_init_block, the CLI)
    (A, "  if (!cfg.config_do_restart || (cfg.config_do_restart && cfg.config_do_DAcycling)) atm_init_coupled_diagnostics(m, s, d);",
     "  atm_init_coupled_diagnostics(m, s, d);",
     RST0("rho_zz pressure_p exner ru rw rtheta_p rho_p") +
     [f"dycore_restart_MW_output60_{f}" for f in "u w rho theta divergence vorticity ke uReconstructZonal uReconstructMeridional".split()]),
    (A, "  if (restart) {", "  if (false) {", RST0("theta_m rho_zz ru rw rtheta_p rho_p exner pressure_p rtheta_base exner_base pressure_base")),
    (A, "      if (r.a->extents() != r.ext)", "      if (false)", ["unit_init_guards"]),
    (RS, "  if (rtrim(xtime) != start_stamp)", "  if (false)", ["err_restart_time", "restart_refusals_MW"]),
    (RS, "  f.Time1 = r.in.state.Time(1);", "", RST0("Time")),
    (X, "    if (cfg.config_do_restart) restart_dycore_fields(rst, f);", "", RST0("Time")),
    (RS, '  return rv && (rv->var_struct == "sfc_input" || rv->var_struct == "diag_physics");',
     '  return rv && rv->var_struct == "sfc_input";', ["dycore_restart0_run_MW"]),
    (RS, '    throw Error(path + ": MPAS reads variable " + n + " from the restart stream; CPPMPAS has no port for it");',
     "    continue;", ["restart_refusals_MW"]),
    (RS, '    if (!f.has_var(n)) throw Error(path + ": variable " + n + " of the restart stream is missing");',
     '    if (false) throw Error(path + ": variable " + n + " of the restart stream is missing");', ["restart_refusals_MW"]),
    (RS, '"pressure_base", "ru", "rw"', '"pressure_base", "rw"', ["dycore_restart0_run_MW"]),
    (X, "    if (!cfg.config_do_restart) check_stream_time(", "    if (false) check_stream_time(", ["restart_refusals_MW"]),
    (X, '    if (cfg.config_do_restart)\n      throw Error("config_do_restart = true: init',
     '    if (false)\n      throw Error("config_do_restart = true: init', ["namelist_MW"]),
    (X, "    if (cfg.config_do_restart) {\n      rst = read_mpas_restart", "    if (false) {\n      rst = read_mpas_restart",
     ["dycore_restart0_run_MW"]),
    (A, "  s.xtime = stamp;                                    // :261-262",
     '  s.xtime = stamp + "x";                              // :261-262', ["dycore_restart0_MW_xtime", "dycore_MW_step0_xtime"]),
    (R, "        in.state.scalars(k, lev, i) = buf[static_cast<std::size_t>(g - 1) * m.nVertLevels + (lev - 1)];",
     "        in.state.scalars(k, lev, i) = buf[static_cast<std::size_t>(g - 1) * m.nVertLevels + (lev - 1)] + 1.0e-12;",
     ["dycore_restart0_MW_qv", "dycore_MW_step0_qv"]),
    ("src/core/namelist_table.inc", r'"config_do_restart", "logical", "false", "false\ntrue")',
     r'"config_do_restart", "logical", "false", "false")', ["dycore_restart0_run_MW", "namelist_MW"]),
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


def test_status(out, t):
    """The status CTest prints on the named test's own line ("Passed", "***Failed", "***Not Run", "***Exception:
    SegFault", "Subprocess aborted***Exception:", ...), or None if the test has no such line."""
    m = re.search(rf"^\s*\d+/\d+ Test +#\d+: {re.escape(t)} \.+ *(.*?)\s+[\d.]+ sec$", out, re.M)
    return m.group(1).strip() if m else None


def killed(status):
    """A test kills the mutant only if it ran and failed. "Not Run" (a fixture setup the mutant broke, a missing
    executable, a disabled test) is a check that did not run: not a kill (CONSTITUTION.md 3.2)."""
    return status is not None and status != "Passed" and "Not Run" not in status


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
            # one CTest run for all the tests the mutant must break: a fixture they share runs once
            r = subprocess.run(["ctest", "--test-dir", str(build), "-R", "^(" + "|".join(map(re.escape, tests)) + ")$",
                                "-FA", "inputs"], capture_output=True, text=True) if tests else None
            for t in tests:
                status = test_status(r.stdout, t)
                failed = r.returncode != 0 and killed(status)
                print(f"{'KILLED' if failed else 'SURVIVED'} {t:22s} by {path}: {new.strip()[:70]}")
                if failed:
                    covered.add(t)
                else:
                    bad.append(f"{t} survived {path}: {new.strip()} (status {status!r})")
        finally:
            p.write_text(orig)
    subprocess.run(["cmake", "--build", str(build), "-j4"], capture_output=True, text=True, check=True)
    for b in bad:
        print("PROBLEM", b)
    print(f"{len(covered)} tests killed by at least one mutant")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
