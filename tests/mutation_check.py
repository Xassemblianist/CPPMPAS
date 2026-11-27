#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Show that every test can fail: apply one deliberate one-line bug at a time to the source, rebuild, run the tests
it must break, require each of them to fail, restore the source. The tests must pass on the unmutated source first
(a test that already fails proves nothing). Run it through the workstation gate.

usage: mutation_check.py <source dir> <build dir> [first | i,j,...]   (index of the first mutant to apply, default 0,
       or a comma-separated list of two or more distinct mutant indices; anything else is refused)

CTest runs one test at a time (-j 1, whatever CTEST_PARALLEL_LEVEL says): tests that share a scratch file must not race in
the check of the unmutated source. A test that times out or is skipped does not kill a mutant. After the last mutant
the source is rebuilt and every test that killed one must pass again, else no kill of the run is trusted.
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
CV = "src/core_atmosphere/diagnostics/mpas_convective_diagnostics.cpp"


def CV0(fields):
    """t = 0 convective diagnostics tests of RG (tests/CMakeLists.txt)"""
    return [f"diag_RG_conv0_{f}" for f in fields.split()]


def CV24(fields, hours="03 06"):
    """3 h / 6 h convective diagnostics tests of RG24 (label long)"""
    return [f"diag_RG24_conv_2026-10-12_{h}.00.00_{f}" for h in hours.split() for f in fields.split()]
RAD, RDRV = "src/physics/module_ra_rrtmg.cpp", "src/physics/mpas_atmphys_radiation.cpp"
RADH = "src/physics/module_ra_rrtmg.hpp"
MC, RF, RD = "src/physics/mpas_atmphys_microphysics.cpp", "src/physics/mp_wsm6_refl.hpp", "src/physics/mp_radar.hpp"


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
    (R, '      in.warnings.push_back("Variable " + n + " not in input file.");',
     '      in.warnings.push_back("Variable " + n + " is not in the input file.");', ["err_missing_scalar"]),
    (R, "  if (vi.dims != want_dims)", "  if (false)", ["err_bad_dims"]),
    (R, "    if (std::find(order.begin(), order.end(), *it) == order.end())", "    if (false)", ["err_unknown_scalar"]),
    (R, "  throw Error(\"variable '\" + name + \"' has no CPPMPAS port (not read or written)\");", "  return r;",
     ["err_unported_field", "unit_io_core"]),
    # a name the Registry does not define: write_history meets field_ref's check first (the message of the CLI,
    # err_unknown_field; without it, field_ref says "has no CPPMPAS port"), registry_var's is tested on its own
    ("src/io/registry.cpp", "  if (!r) throw Error(\"variable '\" + name + \"' is not defined in MPAS-Atmosphere's Registry\");",
     "  if (false) throw Error(\"variable '\" + name + \"' is not defined in MPAS-Atmosphere's Registry\");",
     ["unit_io_core"]),
    (R, "  if (!find_registry_var(name)) throw", "  if (false) throw", ["unit_io_core", "err_unknown_field"]),
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
    (N, "long record) {\n  std::vector<std::size_t> s, c;\n  record_range(name, record, s, c);  // record >= 0 needs a variable whose first dimension is the record dimension\n  if (product(c, path_ + \": write \" + name) != v.size())\n    throw Error(path_ + \": write \" + name + \": size mismatch\");\n  check(nc_put_vara_double(",
     "long record) {\n  std::vector<std::size_t> s, c;\n  record_range(name, -1, s, c); if (record >= 0) { s[0] = record; c[0] = 1; }\n  if (product(c, path_ + \": write \" + name) != v.size())\n    throw Error(path_ + \": write \" + name + \": size mismatch\");\n  check(nc_put_vara_double(",
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
    (W, "    s = s.substr(0, 512);\n", "", ["global_atts"]),
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
    # without config_start_time the clock starts at the Registry default (2010-10-23_00:00:00), which is not the
    # time of the init file: cppmpas run refuses it (check_stream_time, mpas_stream_manager.F:3776-3786) before it
    # writes a step file, so the run fixtures fail and the xtime tests do not run; the xtime values are broken by
    # the value mutants of the writer and of atm_timestep instead
    (NL, "  t(\"config_start_time\", c.config_start_time);", "",
     ["init_step0_MW", "init_step0_SC", "init_step0_JW", "dycore_run_MW", "dycore_run_JW"]),
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
    (A, "    if (n > 25) return kernel_fail", "    if (n > 26) return kernel_fail", ["unit_init_guards"]),
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
    (A, "    if ((nVertLevels % 2 != 0 || nVertLevels <= 2) && std::find(pow_tail.begin(), pow_tail.end(), n) != pow_tail.end()) continue;",
     "    if (false && std::find(pow_tail.begin(), pow_tail.end(), n) != pow_tail.end()) continue;", ["unit_init_guards"]),
    (A, "    if ((nVertLevels % 2 != 0 || nVertLevels <= 2) && std::find(pow_tail.begin(), pow_tail.end(), n) != pow_tail.end()) continue;",
     "    if ((nVertLevels % 2 != 0) && std::find(pow_tail.begin(), pow_tail.end(), n) != pow_tail.end()) continue;", ["unit_init_guards"]),
    # The scalar pow loop for nVertLevels <= 2 (atm_init_coupled_diagnostics) has no mutant: every reference has
    # nVertLevels >= 26, so no test can tell it from the vector loop; the fields it feeds are refused unless
    # --unverified (the two mutants above).
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
    (NL, "  if (kessler || wsm6) v.push_back(\"qr\");",
     "  if (wsm6) v.push_back(\"qr\");", ["unit_namelist", "err_scalars_missing"]),
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
    (A, "    if (n > 25) return kernel_fail", "    if (n >= 25) return kernel_fail", ["unit_init_guards"]),
    (V, "  if (N < 1) return kernel_fail(\"mpas_legs: N < 1\");", "", ["unit_init_guards"]),
    (W, "  if (opt.nprocs < 1 || opt.nprocs > 9) throw", "  if (false) throw", ["unit_io_core"]),
    (W, "    std::remove(path.c_str());", "    (void)path;", ["unit_io_core"]),
    ("src/core/decomp.cpp",
     "  return static_cast<std::size_t>(j - 1) + static_cast<std::size_t>(dim1) * static_cast<std::size_t>(i - 1);",
     "  return static_cast<std::size_t>((j - 1) + dim1 * (i - 1));", ["unit_decomp"]),
    ("src/core/decomp.cpp", "    if (nEdgesOnCell[i - 1] < 0 || nEdgesOnCell[i - 1] > maxEdges)", "    if (false)", ["unit_decomp"]),
    ("src/core/decomp.cpp", "  if (nCellsGlobal < 0) throw", "  if (false) throw", ["unit_decomp"]),
    ("src/core/decomp.cpp", "  if (nCells < 0 || maxCells < 0 || nEdges < 0) throw", "  if (false) throw", ["unit_decomp"]),
    ("src/core/decomp.cpp", "  if (n2 < 0 || static_cast<std::size_t>(n2) > sorted.size())", "  if (false)", ["unit_decomp"]),
    ("src/core/decomp.cpp", "      innerDim > ext(field, 0) || outerDim > ext(field, 1))", "      false)", ["unit_decomp"]),
    ("src/core/decomp.cpp", "  if (innerDim < 0 || outerDim < 0 || indexSpaceDim < 0 || indexSpaceDim > ext(indexToID, 0) ||",
     "  if (innerDim < 0 || outerDim < 0 || indexSpaceDim < 0 ||", ["unit_decomp"]),
    ("src/core/field.hpp", "      if (e != 0 && n > data_.max_size() / static_cast<std::size_t>(e))",
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
    (W, "    return s.substr(0, s.find('\\0'));", "      return s;", ["global_atts"]),
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
    # the harness itself: a test that was not run is not a kill; a test that fails without the mutant is not a kill
    ("tests/mutation_check.py", "    return {t for t in tests if test_status(out, t) != \"Passed\"}",
     "    return {t for t in tests if test_status(out, t) is None}", ["mutation_check_unit"]),
    ("tests/mutation_check.py", "status != \"Passed\" and \"Not Run\" not in status and", "status != \"Passed\" and",
     ["mutation_check_unit"]),
    # --unverified on the CLI (sc_init.nc: nVertLevels = 3) and the exact verified lists
    (X, "    else if (s == \"--unverified\") a.unverified = true;", "    else if (s == \"--unverified\") a.unverified = false;",
     ["init_unverified_SC"]),
    (X, "      if (is_computed && !a.unverified && std::find(verified.begin(), verified.end(), f) == verified.end())",
     "      if (false)", ["err_unverified"]),
    (A, "  const std::vector<std::string> pow_tail = {\"exner\", \"exner_base\", \"pressure_p\", \"pressure_base\"};",
     "  const std::vector<std::string> pow_tail = {\"exner\", \"pressure_p\", \"pressure_base\"};", ["unit_init_guards"]),
    # --- hardening: guards, refusals before the file is created, the separators gfortran
    # fails on, the C locale, the libm barrier of every port
    ("src/core/decomp.cpp", "  const auto ext = [](const FArray<int>& a, int d) { return a.allocated() ? a.extent(d) : 0; };",
     "  const auto ext = [](const FArray<int>& a, int d) { return a.extent(d); };", ["unit_decomp"]),
    ("src/core/decomp.cpp", "      innerDim > ext(field, 0) || outerDim > ext(field, 1))", "      innerDim > ext(field, 0))",
     ["unit_decomp"]),
    ("src/core/decomp.cpp", "  if (nCellsGlobal == INT_MAX) throw", "  if (false) throw", ["unit_decomp"]),
    ("src/core/decomp.cpp", "  if (maxEdges == INT_MAX || nCells == INT_MAX) throw", "  if (maxEdges == INT_MAX) throw",
     ["unit_decomp"]),
    ("src/core/decomp.cpp", "  if (nCells == INT_MAX || maxCells == INT_MAX || nEdges == INT_MAX)",
     "  if (nCells == INT_MAX || maxCells == INT_MAX)", ["unit_decomp"]),
    ("src/core/field.hpp", "      if (e != 0 && n > data_.max_size() / static_cast<std::size_t>(e))",
     "      if (e != 0 && n > std::numeric_limits<std::size_t>::max() / static_cast<std::size_t>(e))", ["unit_decomp"]),
    ("src/core/field.hpp", "    data_.assign(n, fill);  // first:", "    ext_ = extents;\n    data_.assign(n, fill);  // first:",
     ["unit_decomp"]),
    (N, "    if (n > std::numeric_limits<std::size_t>::max() / x)\n      throw", "    if (false)\n      throw", ["unit_io_core"]),
    (A, "    if (ne < 1 || ne > m.maxEdges)", "    if (ne < 1)", ["unit_init_guards"]),
    (A, "    if (ne < 1 || ne > m.maxEdges)", "    if (ne > m.maxEdges)", ["unit_init_guards"]),
    (A, "  if (!m.meshDensity.allocated() || m.meshDensity.extent(0) < m.nCells)\n    throw", "  if (false)\n    throw",
     ["unit_init_guards"]),
    (X, "init_block_verified_fields(in.mesh.nVertLevels, mesh_scaling_unverified(in.mesh, cfg));",
     "init_block_verified_fields(in.mesh.nVertLevels, false);", ["err_scaled_mesh"]),
    (NL, "  return strtod_l(s.c_str(), nullptr, c);", "  return (void)c, std::strtod(s.c_str(), nullptr);", ["unit_namelist_locale"]),
    (NL, "    if (commas > 1) throw", "    if (false) throw", ["unit_namelist"]),
    (NL, "    if (v.quoted && i < s.size() && !std::isspace", "    if (false && i < s.size() && !std::isspace", ["unit_namelist"]),
    (NL, "    out[key] = match;", "    out[key] = v;", ["unit_namelist"]),
    (NL, "  } else {\n    throw Error(\"active scalars: config_physics_suite", "  } else if (false) {\n    throw Error(\"active scalars: config_physics_suite", ["unit_namelist"]),
    (NL, "  if (c.config_les_model != \"none\") throw", "  if (false) throw", ["unit_namelist"]),
    (NL, "  if (cu != \"off\" && cu != \"cu_ntiedtke\") throw", "  if (false) throw", ["unit_namelist"]),
    (NL, "  if (pbl != \"off\" && pbl != \"bl_ysu\") throw", "  if (false) throw", ["unit_namelist"]),
    ("tools/gen_namelist_table.py", "            if rec in seen:", "            if False:", ["namelist_gen_errors"]),
    (W, "    if (!names.insert(n).second) throw", "    if (false) throw", ["unit_io_core"]),
    (W, "    if (ref.c && ref.c->size() > kStrLen) throw", "    if (false) throw", ["unit_io_core"]),
    (W, "      if (d != \"Time\") (void)file_dim(m, d);", "      (void)d;", ["unit_io_core"]),
    (W, "    f.put_att_text(\"\", \"git_version\", trimmed(opt.git_version));",
     "    f.put_att_text(\"\", \"git_version\", opt.git_version);", ["unit_io_core"]),
    (W, "    while (!s.empty() && s.back() == ' ') s.pop_back();\n    return s.substr(0, s.find('\\0'));",
     "      return s.substr(0, s.find('\\0'));", ["global_atts"]),
    # an empty parent_id not written at all: only a check of NetCDF's status finds the attribute missing
    (W, "    f.put_att_text(\"\", \"parent_id\", trimmed(m.parent_id));",
     "    if (!trimmed(m.parent_id).empty()) f.put_att_text(\"\", \"parent_id\", trimmed(m.parent_id));", ["global_atts"]),
    (R, "rtrim(file_id).empty() ? parent_id :", "false ? parent_id :", ["global_atts"]),
    (R, "(rtrim(file_id) + \"\\n\" + rtrim(parent_id)).substr(0, kStrKIND);", "(rtrim(file_id) + \"\\n\" + rtrim(parent_id));",
     ["unit_io_core"]),
    (L, "CPPMPAS_FORCEINLINE double sin(double x) { return ::cppwrf::glibc::sin(arg(x)); }",
     "CPPMPAS_FORCEINLINE double sin(double x) { return ::cppwrf::glibc::sin(x); }", ["unit_libm"]),
    (L, "CPPMPAS_FORCEINLINE double cos(double x) { return ::cppwrf::glibc::cos(arg(x)); }",
     "CPPMPAS_FORCEINLINE double cos(double x) { return ::cppwrf::glibc::cos(x); }", ["unit_libm"]),
    # removing arg() from asin, acos, atan2, pow or mvec_pow (either argument) changes no bit under g++ 14.2 -O3,
    # even with -x written at the call (unit_libm calls every port so): g++ does not fold the negation there, so those
    # mutants are equivalent today; unit_libm would fail on a compiler that folds it
    (L, "  x = arg(x);\n  constexpr double hp0", "  constexpr double hp0", ["unit_libm"]),
    # --- the harness: what kills, the selection, the baseline in one process at a time, the restored
    # source checked again (unit tests of main() with a simulated cmake / ctest)
    ("tests/mutation_check.py", " and \"Skipped\" not in status\n", "\n", ["mutation_check_unit"]),
    ("tests/mutation_check.py", "            and \"Timeout\" not in status)", "            )", ["mutation_check_unit"]),
    ("tests/mutation_check.py", "\n    return rc != 0 and killed(status) and t not in broken",
     "\n    return rc != 0 and killed(status)", ["mutation_check_unit"]),
    ("tests/mutation_check.py", "\n    return rc != 0 and killed(status) and t not in broken",
     "\n    return killed(status) and t not in broken", ["mutation_check_unit"]),
    ("tests/mutation_check.py", "    if not re.fullmatch(r\"[0-9]+(,[0-9]+)+\", sel):", "    if not re.match(r\"[0-9]+(,[0-9]+)+\", sel):",
     ["mutation_check_unit"]),
    ("tests/mutation_check.py", "\n    if len(set(idx)) != len(idx):", "\n    if False:", ["mutation_check_unit"]),
    ("tests/mutation_check.py", "        if i >= n:\n", "        if i > n:\n", ["mutation_check_unit"]),
    ("tests/mutation_check.py", "\n        if int(sel) >= n:", "\n        if int(sel) > n:", ["mutation_check_unit"]),
    ("tests/mutation_check.py", "[\"ctest\", \"--test-dir\", str(build), \"-j\", \"1\", \"-R\"", "[\"ctest\", \"--test-dir\", str(build), \"-R\"",
     ["mutation_check_unit"]),
    ("tests/mutation_check.py", "    return \"^(\" + \"|\".join(re.escape(t) for t in tests) + \")$\"",
     "    return \"^(\" + \"|\".join(tests) + \")$\"", ["mutation_check_unit"]),
    ("tests/mutation_check.py", "\n        for t in sorted(not_passing(r.stdout, covered)):", "\n        for t in []:",
     ["mutation_check_unit"]),
]


# ---------------------------------------------------------------------------------------------------------------
# the dycore (src/core_atmosphere/atm_dycore.cpp, atm_dissipation.cpp) against the MW step files, its checks and the
# run command
DY, DS = "src/core_atmosphere/atm_dycore.cpp", "src/core_atmosphere/atm_dissipation.cpp"


def ST(steps, fields="u w rho_zz theta_m rho theta pressure_p"):
    """dycore_MW_step<n>_<field> tests (tests/CMakeLists.txt)"""
    return [f"dycore_MW_step{n}_{f}" for n in steps for f in fields.split()]


def JST(steps, fields="u w rho_zz theta_m rho theta pressure_p"):
    """dycore_JW_step<n>_<field> tests (tests/CMakeLists.txt)"""
    return [f"dycore_JW_step{n}_{f}" for n in steps for f in fields.split()]


S1 = ST([1])
ALLDYN = ST(range(1, 11)) + ["dycore_MW_mpas_compare"]
LONG_R = "u w rho_zz theta_m pressure_p rho theta exner ru rw rtheta_p rho_p".split()
LONG_O = "u w pressure rho theta divergence vorticity ke uReconstructZonal uReconstructMeridional".split()
LONG = ([f"dycore_MW_restart_{f}" for f in LONG_R] + [f"dycore_MW_output30_{f}" for f in LONG_O] +
        [f"dycore_MW_output60_{f}" for f in LONG_O])
J1 = JST([1])
JALLDYN = JST(range(1, 11)) + ["dycore_JW_mpas_compare"]
JLONG = [f"dycore_JW_restart_{f}" for f in LONG_R] + [f"dycore_JW_output24h_{f}" for f in LONG_O]
MUTATIONS += [
    # atm_rk_integration_setup, moist coefficients, physics_get_tend, atm_compute_vert_imp_coefs
    (DY, "      f.u2(k, iEdge) = s.u(k, iEdge);", "      f.u2(k, iEdge) = 0.999 * s.u(k, iEdge);", S1),
    (DY, "// :3165-3171\n    for (int k = 1; k <= nVertLevels + 1; ++k) {\n      f.rw_save(k, iCell) = d.rw(k, iCell);",
     "// :3165-3171\n    for (int k = 1; k <= nVertLevels + 1; ++k) {\n      f.rw_save(k, iCell) = 0.0;", S1),
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
    # the clamp of TI:2799-2800 acts on no reference (MW's qv and SC's mono-limited scalars are >= +0 there): the
    # call deleted is exact on every step file; unit_dycore tests the routine itself
    (DY, "    if (scalars_2[i] < 0.0) scalars_2[i] = 0.0;", "    ;", ["unit_dycore"]),
    (DY, "    if (scalars_2[i] < 0.0) scalars_2[i] = 0.0;", "    if (scalars_2[i] <= 0.0) scalars_2[i] = 0.0;",
     ["unit_dycore"]),
    (DY, "  parallel_for_n(scalars.size(), [=] CPPMPAS_LAMBDA(std::size_t i) {",
     "  parallel_for_n(scalars.size() - 3, [=] CPPMPAS_LAMBDA(std::size_t i) {", ["unit_dycore"]),
    (DY, "    if (!moist) continue;", "", ["unit_dycore"]),
    # output diagnostics, Time, xtime, the step files
    (DY, "      s.rho(k, iCell) = s.rho_zz(k, iCell) * m.zz(k, iCell);", "      s.rho(k, iCell) = s.rho_zz(k, iCell);",
     ST(range(11), "rho") + JST(range(11), "rho")),
    (DY, "      s.theta(k, iCell) = s.theta_m(k, iCell) / (1. + rvord * s.scalars(index_qv, k, iCell));",
     "      s.theta(k, iCell) = s.theta_m(k, iCell) / (1.5 + rvord * s.scalars(index_qv, k, iCell));", ST(range(11), "theta") + JST(range(11), "theta")),
    (DY, "  f.Time2 = s.Time(1) + dt;", "  f.Time2 = s.Time(1) + dt + 1.0;", ST(range(1, 11), "Time") + ["dycore_MW_restart_Time"] +
     JST(range(1, 11), "Time") + ["dycore_JW_restart_Time"]),
    # Time read from the input stream (review F3); the copy into time level 2 at allocation (mpas_atm_core.F:213-218)
    # has no mutant: atm_timestep overwrites Time2 before anything reads it
    (R, '"theta_base", "surface_pressure", "Time"};', '"theta_base", "surface_pressure"};', ["dycore_run_MW", "unit_dycore"]),
    (DY, "  std::swap(s.Time(1), f.Time2);", "", ST(range(1, 11), "Time") + JST(range(1, 11), "Time") + ["unit_dycore"]),
    (X, "      atm_timestep(dcfg, in.mesh, in.state, in.diag, f, dcfg.config_dt, add_seconds(start, n * dt_s));",
     "      atm_timestep(dcfg, in.mesh, in.state, in.diag, f, dcfg.config_dt, add_seconds(start, (n + 1) * dt_s));",
     ST([1], "xtime u w") + ["dycore_restart_MW_output60_xtime"]),
    (X, "    write_step(0);", "", ST([0], "xtime Time u w rho_zz theta_m rho theta pressure_p qv") +
     JST([0], "xtime Time u w rho_zz theta_m rho theta pressure_p qv")),
    (X, "    for (const auto& sn : in.state.scalar_names) fields.push_back(sn);", "",
     ST(range(11), "qv") + JST(range(11), "qv")),
    ("src/io/mpas_writer.cpp", '          f.put_att_text(n, k, time_units);', "          f.put_att_text(n, k, v);",
     ST(range(11), "Time") + JST(range(11), "Time")),
    # the init pow loop for nVertLevels <= 2: no reference has such a mesh, so its fields need
    # --unverified; the scalar loop itself (atm_init.cpp) has no mutant for the same reason
    (A, "    if ((nVertLevels % 2 != 0 || nVertLevels <= 2) &&", "    if ((nVertLevels % 2 != 0) &&", ["unit_init_guards"]),
    # mutants left out earlier as exact on MW that JW kills: the level index of etp / ewp (JW's top at
    # 45 km crosses the epssm transition, so etp / ewp vary with k), sign(1, ru) against a constant +1 (JW's flow
    # reverses)
    (DY, "      a_tri(k, iCell) = a_tri(k, iCell) * etp(k - 1) * ewp(k - 1);",
     "      a_tri(k, iCell) = a_tri(k, iCell) * etp(k) * ewp(k - 1);", J1),
    (DY, "      c_tri[k] = c_tri[k] * etp(k) * ewp(k + 1);", "      c_tri[k] = c_tri[k] * etp(k) * ewp(k);", J1),
    (DY, "m.adv_coefs(j, iEdge) + sign1(ru(k, iEdge)) * m.adv_coefs_3rd(j, iEdge);",
     "m.adv_coefs(j, iEdge) + 1.0 * m.adv_coefs_3rd(j, iEdge);", J1),
    (DY, "m.adv_coefs(j, iEdge) + sign1(ru_edge_w[k]) * m.adv_coefs_3rd(j, iEdge);",
     "m.adv_coefs(j, iEdge) + 1.0 * m.adv_coefs_3rd(j, iEdge);", J1),
    # Exact on MW and JW for the reasons given, so not in the list (the first three were run on JW, the rk_timestep,
    # ruAvg, rthdynten and rho_zz_old_split ones on MW):
    #  * the reference_u lanes all sincos ("vec = false") or the reference_u term deleted: the term is +-0 on both
    #    (fEdge = 0 on MW, u_init = v_init = 0 on JW); all 71 JW tests pass with either mutant;
    #  * the smagorinsky_2d bound (0.01*len**2)*invDt at 1/100 of its value: it never binds on JW;
    #  * rk_timestep[1] / rk_timestep[2]: atm_recover_large_step_variables reads dt only at rk_step 3 (TI:4437) and
    #    the other readers (advance_scalars, LBCs) are not on the ported path (the earlier reason, rt_diabatic_tend = 0,
    #    was wrong);
    #  * ruAvg / wwAvg (their invNs averaging and the _split copies): read only by advance_scalars and split
    #    transport, not ported here; rthdynten: read only by convection schemes; rho_zz_old_split: copied into
    #    time level 1 at substep_finish, which nothing reads before the level swap and the next setup overwrite it;
    #  * qtot sums and the rb*qtot terms: qv = 0 on MW and JW (SC exercises them);
    #  * the copy of Time into time level 2 at allocation: atm_timestep overwrites it before it is read.
    # JW: smagorinsky_2d (DM:119-204) and its call (TI:6346-6352)
    (DS, "        dudx[k] = dudx[k] + c2 * u(k, e) - cs * v(k, e);", "        dudx[k] = dudx[k] + c2 * u(k, e) + cs * v(k, e);",
     JALLDYN + JLONG),
    (DS, "        dudy[k] = dudy[k] + cs * u(k, e) - s2 * v(k, e);", "        dudy[k] = dudy[k] + cs * u(k, e) - cs * v(k, e);", J1),
    (DS, "        dvdx[k] = dvdx[k] + cs * u(k, e) + c2 * v(k, e);", "        dvdx[k] = dvdx[k] + cs * u(k, e) + c2 * u(k, e);", J1),
    (DS, "        dvdy[k] = dvdy[k] + s2 * u(k, e) + cs * v(k, e);", "        dvdy[k] = dvdy[k] + s2 * u(k, e) - cs * v(k, e);", J1),
    (DS, "      const real d_12 = dudy[k] + dvdx[k];", "      const real d_12 = dudy[k] - dvdx[k];", J1),
    (DS, "      const real d_22 = 2 * dvdy[k];", "      const real d_22 = 2 * dudy[k];", J1),
    (DS, "kfac * std::sqrt(0.25 * ((d_11 - d_22) * (d_11 - d_22)) + d_12 * d_12);",
     "kfac * std::sqrt(0.5 * ((d_11 - d_22) * (d_11 - d_22)) + d_12 * d_12);", J1),
    (DS, "  const real kfac = cl * cl;", "  const real kfac = cl;", J1),
    (DS, "      const real c2 = m.deformation_coef_c2(iEdge, iCell), s2 = m.deformation_coef_s2(iEdge, iCell);",
     "      const real c2 = m.deformation_coef_s2(iEdge, iCell), s2 = m.deformation_coef_s2(iEdge, iCell);", J1),
    (DS, "      dudx[k] = 0.0;\n", "\n", J1),
    (DS, "  h_mom_eddy_visc4 = config_visc4_2dsmag * (config_len_disp * config_len_disp * config_len_disp);",
     "  h_mom_eddy_visc4 = config_visc4_2dsmag * (config_len_disp * config_len_disp);", J1),
    (DS, "  h_theta_eddy_visc4 = h_mom_eddy_visc4;", "  h_theta_eddy_visc4 = 0.0;", JST([1], "theta_m rho_zz theta")),
    (DY, "      smagorinsky_2d(m, f.eddy_visc_horz, f.u2, d.v, cfg.config_smagorinsky_coef,",
     "      smagorinsky_2d(m, f.eddy_visc_horz, f.u2, f.u2, cfg.config_smagorinsky_coef,", J1),
    # JW: the 4th-order filters (DM:728-811, 1073-1100, 1306-1332)
    (DS, "        const real edge_sign = m.invAreaTriangle(iVertex) * m.dcEdge(iEdge) * m.edgesOnVertex_sign(i, iVertex);",
     "        const real edge_sign = m.invAreaTriangle(iVertex) * m.dvEdge(iEdge) * m.edgesOnVertex_sign(i, iVertex);", J1),
    (DS, "      for (int k = 1; k <= nVertLevels; ++k) delsq_vorticity(k, iVertex) = 0.0;", "", JST([2])),
    (DS, "      for (int k = 1; k <= nVertLevels; ++k) delsq_divergence(k, iCell) = 0.0;", "", JST([2])),
    (DS, "        const real edge_sign = r * m.dvEdge(iEdge) * m.edgesOnCell_sign(i, iCell);",
     "        const real edge_sign = r * m.dcEdge(iEdge) * m.edgesOnCell_sign(i, iCell);", J1),
    (DS, "      const real r_dc = u_mix_scale * a.config_del4u_div_factor * m.invDcEdge(iEdge);",
     "      const real r_dc = u_mix_scale * m.invDcEdge(iEdge);", J1),
    (DS, "      const real r_dv = u_mix_scale * std::min(m.invDvEdge(iEdge), 4 * m.invDcEdge(iEdge));",
     "      const real r_dv = u_mix_scale * std::min(m.invDvEdge(iEdge), 1 * m.invDcEdge(iEdge));", J1),
    (DS, "        tend_u_euler(k, iEdge) = tend_u_euler(k, iEdge) - u_diffusion;  // DM:804",
     "        tend_u_euler(k, iEdge) = tend_u_euler(k, iEdge) + u_diffusion;  // DM:804", J1),
    (DS, "          tend_w_euler(k, iCell) = tend_w_euler(k, iCell) - edge_sign * (delsq_w(k, cell2) - delsq_w(k, cell1));",
     "          tend_w_euler(k, iCell) = tend_w_euler(k, iCell) + edge_sign * (delsq_w(k, cell2) - delsq_w(k, cell1));",
     JST([1], "w")),
    (DS, "      const real r_areaCell = a.h_mom_eddy_visc4 * m.invAreaCell(iCell);",
     "      const real r_areaCell = a.h_mom_eddy_visc4 * m.invAreaCell(iCell) * 0.5;", JST([1], "w")),
    (DS, "              tend_theta_euler(k, iCell) - edge_sign * (delsq_theta(k, cell2) - delsq_theta(k, cell1));",
     "              tend_theta_euler(k, iCell) - edge_sign * (delsq_theta(k, cell2) - delsq_theta(k, cell2));",
     JST([1], "theta_m theta")),
    (DS, "      const real r_areaCell = a.h_theta_eddy_visc4 * a.prandtl_inv * m.invAreaCell(iCell);",
     "      const real r_areaCell = a.h_theta_eddy_visc4 * a.prandtl_inv * m.invAreaCell(iCell) * 0.5;",
     JST([1], "theta_m theta")),
    (DY, "  da.config_del4u_div_factor = cfg.config_del4u_div_factor;", "  da.config_del4u_div_factor = 1.0;", J1),
    # step 0 by value with the file names intact (review F12): every real field written as 1.0 (w is 0 at t = 0),
    # character fields padded with 'x'
    ("src/io/mpas_writer.cpp", "    std::vector<double> local(ref.r->data(), ref.r->data() + ref.r->size());",
     "    std::vector<double> local(ref.r->size(), 1.0);",
     ST([0], "Time u w rho_zz theta_m rho theta pressure_p") + JST([0], "Time u w rho_zz theta_m rho theta pressure_p")),
    ("src/io/mpas_writer.cpp", "    v.resize(kStrLen, ' ');", "    v.resize(kStrLen, 'x');",
     ["dycore_MW_step0_xtime", "dycore_JW_step0_xtime"]),
    # the step files' global attributes (review F14): tests/dycore_compare.py checks them file by file
    ("src/io/mpas_writer.cpp", '    f.put_att_text("", "version", "8.4.2");', '    f.put_att_text("", "version", "8.4.1");',
     ["dycore_MW_mpas_compare", "dycore_JW_mpas_compare"]),
    ("src/io/mpas_writer.cpp", '    f.put_att_double("", "sphere_radius", m.sphere_radius);',
     '    f.put_att_double("", "sphere_radius", m.sphere_radius * (1 + 1e-15));', ["dycore_MW_mpas_compare", "dycore_JW_mpas_compare"]),
    # output diagnostics (review F10): pressure (mpas_atm_core.F:961); computing them only on output steps has no
    # mutant (theta, rho and pressure are read by no routine of the step)
    (DY, "      d.pressure(k, iCell) = d.pressure_base(k, iCell) + d.pressure_p(k, iCell);",
     "      d.pressure(k, iCell) = d.pressure_p(k, iCell);", ["dycore_MW_output30_pressure", "dycore_MW_output60_pressure", "dycore_JW_output24h_pressure"]),
    # time stamps (review F8, F18); the blank padding of xtime itself has no file-level mutant: the writer pads
    # character fields with blanks too (mpas_writer.cpp write_values), unit_dycore compares the string
    (DY, "    form = want == '0' ? (ts[i] >= '0' && ts[i] <= '9') : ts[i] == want;",
     "    form = want == '0' ? true : ts[i] == want;", ["unit_dycore"]),
    (DY, "  if (seconds > 9999LL * 366 * 86400) throw", "  if (false) throw", ["unit_dycore"]),
    (DY, '  if (Y > 9999) throw Error("add_seconds:', '  if (false) throw Error("add_seconds:', ["unit_dycore"]),
    (DY, '  if (xtime_new.size() > 64) throw', '  if (false) throw', ["unit_dycore"]),
    (DY, "  f.xtime2.resize(64, ' ');", "", ["unit_dycore"]),
    # the content of xtime with its file name unchanged (review F12): the files' xtime values must differ
    (DY, "  f.xtime2.resize(64, ' ');", "  f.xtime2.resize(64, 'x');",
     ST(range(1, 11), "xtime") + JST(range(1, 11), "xtime") + ["dycore_MW_restart_xtime", "dycore_MW_output30_xtime",
      "dycore_MW_output60_xtime", "dycore_JW_restart_xtime", "dycore_JW_output24h_xtime", "unit_dycore"]),
    ("src/io/mpas_writer.cpp", '    if (t.empty()) throw Error("write Time: initial_time was not read");',
     '    if (false) throw Error("write Time: initial_time was not read");', ["unit_dycore"]),
    # refusals and array checks of the dycore review (F4, F5, F6, F17); unit_dycore kills each
    (DY, "  if (!(c.config_mpas_cam_coef == 0.0)) no(", "  if (c.config_mpas_cam_coef > 0.0) no(", ["unit_dycore"]),
    (DY, "  if (!(c.config_h_mom_eddy_visc4 == 0.0)) no(", "  if (c.config_h_mom_eddy_visc4 > 0.0) no(", ["unit_dycore"]),
    (DY, "  if (!(c.config_h_theta_eddy_visc4 == 0.0)) no(", "  if (false) no(", ["unit_dycore"]),
    (DY, "  if (!(c.config_apvm_upwinding > 0.0)) no(", "  if (false) no(", ["unit_dycore"]),
    (DY, "    if (!std::isfinite(x)) no(", "    if (std::isnan(x)) no(", ["unit_dycore"]),
    (DY, "  if (!off(c.config_lsm_scheme)) no(", "  if (false) no(", ["unit_dycore"]),
    (DY, "  if (!off(c.config_sfclayer_scheme)) no(", "  if (false) no(", ["unit_dycore"]),
    (DY, "  if (!off(c.config_gwdo_scheme)) no(", "  if (false) no(", ["unit_dycore"]),
    (DY, "  if (c.config_oml1d) no(", "  if (false) no(", ["unit_dycore"]),
    (DY, '  if (c.config_calendar_type != "gregorian") no(', "  if (false) no(", ["unit_dycore"]),
    (DY, " || !off(c.config_radt_cld_scheme))", ")", ["unit_dycore"]),
    (DY, "      if (m.meshDensity(iCell) != 1.0) no(", "      if (false) no(", ["unit_dycore"]),
    (DY, "  c.config_h_ScaleWithMesh = init_cfg.config_h_ScaleWithMesh;", "  c.config_h_ScaleWithMesh = true;", ["unit_dycore"]),
    (DY, "  if (s.u_init.rank() != 1 || s.u_init.extent(0) < m.nVertLevels ||", "  if (false &&", ["unit_dycore"]),
    (DY, "    if (m.bdyMaskCell(iCell) != 0 || m.specZoneMaskCell(iCell) != 0.0) bad(",
     "    if (m.bdyMaskCell(iCell) != 0) bad(", ["unit_dycore"]),
    (DY, "    if (m.bdyMaskEdge(iEdge) != 0 || m.specZoneMaskEdge(iEdge) != 0.0) bad(",
     "    if (m.specZoneMaskEdge(iEdge) != 0.0) bad(", ["unit_dycore"]),
    (DY, "  for (const auto* a : {&f.ru_save, &f.cqu, &f.ru_p,", "  for (const auto* a : {&f.ru_save, &f.ru_p,", ["unit_dycore"]),
    (DY, "has(m.zgrid, {nL + 1, nCs}, false)", "has(m.zgrid, {nL, nCs}, false)", ["unit_dycore"]),
    (DY, "    if (m.nAdvCellsForEdge(iEdge) < 0 || m.nAdvCellsForEdge(iEdge) > m.advCellsForEdge.extent(0) ||",
     "    if (m.nAdvCellsForEdge(iEdge) < 0 ||", ["unit_dycore"]),
    (DY, "    if (m.nEdgesOnCell(iCell) < 0 || m.nEdgesOnCell(iCell) > mE) bad(", "    if (false) bad(", ["unit_dycore"]),
    (DY, "  if (nL < 4) bad(", "  if (false) bad(", ["unit_dycore"]),
    (DY, " || f.moist_end > nsf ||", " ||", ["unit_dycore"]),
    (DY, "  check_arrays(m, &s, &d, f);\n  const int nVertLevels = m.nVertLevels, num_scalars",
     "  const int nVertLevels = m.nVertLevels, num_scalars", ["unit_dycore"]),
    (DY, "void atm_compute_moist_coefficients(const Mesh& m, DycoreFields& f) {\n  check_arrays(m, nullptr, nullptr, f);",
     "void atm_compute_moist_coefficients(const Mesh& m, DycoreFields& f) {", ["unit_dycore"]),
    (DY, "void physics_get_tend(const Mesh& m, const State& s, DycoreFields& f) {\n  check_arrays(m, &s, nullptr, f);",
     "void physics_get_tend(const Mesh& m, const State& s, DycoreFields& f) {", ["unit_dycore"]),
    (DY, "DycoreFields& f, real dts) {\n  check_arrays(m, &s, &d, f);", "DycoreFields& f, real dts) {", ["unit_dycore"]),
    (DY, "  check_dycore_supported(cfg, m, s);\n  check_arrays(m, &s, &d, f);\n  if (rk_step < 1",
     "  check_dycore_supported(cfg, m, s);\n  if (rk_step < 1", ["unit_dycore"]),
    (DY, "void atm_set_smlstep_pert_variables(const Mesh& m, DycoreFields& f) {\n  check_arrays(m, nullptr, nullptr, f);",
     "void atm_set_smlstep_pert_variables(const Mesh& m, DycoreFields& f) {", ["unit_dycore"]),
    (DY, "  check_arrays(m, &s, &d, f);\n  if (dynamics_split < 1 ||", "  if (dynamics_split < 1 ||", ["unit_dycore"]),
    (DY, '  if (rk_step < 1 || rk_step > 3) throw Error("atm_compute_dyn_tend:', '  if (false) throw Error("atm_compute_dyn_tend:',
     ["unit_dycore"]),
    (DY, "  if (rk_step < 1 || rk_step > 3 || ns < 1) throw", "  if (rk_step < 1 || rk_step > 3) throw", ["unit_dycore"]),
    (DY, '  if (small_step < 1) throw', '  if (false) throw', ["unit_dycore"]),
    (DY, "  if (dynamics_split < 1 || dynamics_substep < 1 || dynamics_substep > dynamics_split)\n", "  if (false)\n", ["unit_dycore"]),
    (DY, "      m.nVertices > kMax || m.nVertLevels > kMax)\n", "      false)\n", ["unit_dycore"]),
    (DY, "  if (s.scalars.rank() != 3 || s.scalars.extent(0) < 1 ||", "  if (false &&", ["unit_dycore"]),
    # not mutated: the check_arrays calls of atm_advance_acoustic_step, atm_divergence_damping_3d,
    # atm_recover_large_step_variables and atm_compute_output_diagnostics' check (deleting one makes the unit test
    # write out of bounds, ru_p / theta); atm_srk3's own check_arrays (each routine it calls checks the same arrays
    # first: deleting it changes nothing)
    # the configuration checks, the CLI
    (DY, '  if (c.config_horiz_mixing != "2d_fixed" && c.config_horiz_mixing != "2d_smagorinsky")',
     '  if (c.config_horiz_mixing != "2d_fixed")', ["dycore_run_JW", "dycore_run_long_JW", "unit_dycore"]),
    (DY, '  if (c.config_horiz_mixing != "2d_fixed" && c.config_horiz_mixing != "2d_smagorinsky")',
     '  if (c.config_horiz_mixing != "2d_smagorinsky")',
     ["dycore_run_MW", "dycore_run_long_MW", "dycore_restart0_run_MW", "dycore_restart_run_long_MW", "unit_dycore"]),
    (DY, '  if (c.config_horiz_mixing != "2d_fixed" && c.config_horiz_mixing != "2d_smagorinsky")',
     '  if (c.config_horiz_mixing.rfind("2d_", 0) != 0)', ["unit_dycore"]),
    (DY, "  if (c.config_scalar_advection) check_advance_scalars_supported(c, m);",
     "  if (false) check_advance_scalars_supported(c, m);", ["unit_dycore"]),
    (DY, "      if (m.fEdge(iEdge) != 0.0) no(", "      if (false) no(", ["unit_dycore"]),
    (DY, "  if (!uv_init_zero)\n", "  if (true)\n", ["dycore_run_JW", "dycore_run_long_JW", "unit_dycore"]),
    (DY, "    if (s.u_init(k) != 0.0 || s.v_init(k) != 0.0) uv_init_zero = false;",
     "    if (s.u_init(k) != 0.0) uv_init_zero = false;", ["unit_dycore"]),
    (DY, "  for (int k = 1; k <= m.nVertLevels; ++k)\n    if (s.u_init(k)", "  for (int k = 1; k < m.nVertLevels; ++k)\n    if (s.u_init(k)",
     ["unit_dycore"]),
    (DY, "  if (m.nVertLevels < 4) no(", "  if (false) no(", ["err_dycore_odd_levels", "unit_dycore"]),
    (DY, "  if (c.config_time_integration_order != 2)", "  if (c.config_time_integration_order > 3)", ["unit_dycore"]),
    (DY, "  if (!(c.config_v_mom_eddy_visc2 >= 0.0)) no(", "  if (false) no(", ["unit_dycore"]),
    (DY, "  if (!(c.config_v_mom_eddy_visc2 >= 0.0)) no(", "  if (!(c.config_v_mom_eddy_visc2 > 0.0)) no(",
     ["dycore_run_JW", "dycore_run_long_JW", "unit_dycore"]),
    # (the NaN-blind form "< 0.0" is exact now: a NaN viscosity is refused earlier as a real that is not finite)
    (DY, "  if (!(c.config_v_theta_eddy_visc2 >= 0.0)) no(", "  if (false) no(", ["unit_dycore"]),
    (DY, '  c.config_smagorinsky_coef = r("config_smagorinsky_coef");', '  c.config_smagorinsky_coef = r("config_visc4_2dsmag");',
     J1 + ["unit_dycore"]),
    (DY, '  c.config_visc4_2dsmag = r("config_visc4_2dsmag");', '  c.config_visc4_2dsmag = r("config_smagorinsky_coef");',
     J1 + ["unit_dycore"]),
    (DY, "  if (!off(c.config_microp_scheme) && !kessler) no(", "  if (false) no(", ["unit_dycore"]),
    (DY, "  if (!(c.config_dt > 0.0) || c.config_dt != std::floor(c.config_dt) || c.config_dt > 86400.0)",
     "  if (!(c.config_dt > 0.0) || c.config_dt > 86400.0)", ["unit_dycore"]),
    (DY, "  c.config_len_disp = init_cfg.config_len_disp;", '  c.config_len_disp = r("config_len_disp");', ["unit_dycore"]),
    (DY, "  auto leap = [](int y) { return (y % 4 == 0 && y % 100 != 0) || y % 400 == 0; };",
     "  auto leap = [](int y) { return y % 4 == 0; };", ["unit_dycore"]),
    (DY, "    if (++D > mdays(Y, M)) {", "    if (++D >= mdays(Y, M)) {", ["unit_dycore"]),
    (X, "  if (used != v.size() || n < 0 || n > 1000000) throw Error(opt +", "  if (n < 0 || n > 1000000) throw Error(opt +",
     ["err_dycore_steps"]),
    (DY, "  mpas_reconstruct_2d(m, f.u2, d.uReconstructX,", "  if (false) mpas_reconstruct_2d(m, f.u2, d.uReconstructX,",
     ["dycore_MW_output30_uReconstructZonal", "dycore_MW_output30_uReconstructMeridional", "dycore_MW_output60_uReconstructZonal"]),
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
    (RS, '  if (!NcFile::open_read(path).has_var("Time")) throw', '  if (false) throw', ["restart_refusals_MW"]),
    # the restart's Time is State::Time(1) as read (the copy into f.Time2 is overwritten before it is read: no mutant)
    (RS, '  read_field_ref(f, in, "Time");', '  read_field_ref(f, in, "Time"); in.state.Time(1) = 0.0;', RST0("Time")),
    # (deleting the restart_dycore_fields call is exact since Time lives in State::Time: ru_p and rw_p are overwritten
    # at small step 1, rt_diabatic_tend is 0 on MW and rthdynten has no reader on the ported path)
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
    # the CPU build has none of the GPU's capacity limits (core/backend.hpp, operators/vector_operations.cpp)
    ("src/core/backend.cpp", "#if !defined(CPPMPAS_USE_CUDA)\n  (void)nVertLevels;", "#if 0\n  (void)nVertLevels;",
     ["unit_backend"]),
    ("src/core/backend.hpp", "    if (n_ > N) {\n      big_.assign", "    if (n_ > N) {\n      n_ = N;\n      big_.assign",
     ["unit_backend"]),
    (V, "#if defined(CPPMPAS_USE_CUDA)\n  if (N > kMaxRbfPoints + 2)", "#if 1\n  if (N > kMaxRbfPoints + 2)",
     ["unit_backend"]),
    (V, "#if defined(CPPMPAS_USE_CUDA)\n  if (pointCount > kMaxRbfPoints)", "#if 1\n  if (pointCount > kMaxRbfPoints)",
     ["unit_backend"]),
]

# ---------------------------------------------------------------------------------------------------------------
# the supercell: scalar transport (atm_advance_scalars.cpp), perturbation mixing (atm_dissipation_perturbation.cpp),
# Kessler microphysics (atm_microphysics.cpp, src/physics/mp_kessler.hpp) and their wiring, against the SC step files
AS, DP = "src/core_atmosphere/atm_advance_scalars.cpp", "src/core_atmosphere/atm_dissipation_perturbation.cpp"
MP, KS = "src/core_atmosphere/atm_microphysics.cpp", "src/physics/mp_kessler.hpp"
SCF = "u w rho_zz theta_m rho theta pressure_p qv qc qr"


def SCT(steps, fields=SCF):
    """dycore_SC_step<n>_<field> tests (tests/CMakeLists.txt)"""
    return [f"dycore_SC_step{n}_{f}" for n in steps for f in fields.split()]


# qc and qr are 0 in every cell of the SC reference's first 30 s (no cloud yet): the transport and Kessler mutants
# are required to break qv (and theta_m / pressure_p where named); qc and qr have their own mutants below
SCQ = lambda steps: SCT(steps, "qv")
# the Kessler kernel's condensation, rain and sedimentation act only once cloud forms: their mutants are exact on the
# ten steps (every column subsaturated: product = max(prod, -qc) = 0, qrprod = ern = 0) and must break the 30-min
# state (label long)
# the limiter acts where the high-order update would leave the local bounds, which the smooth qv of the first 30 s
# does not reach: these mutants must break the 30-min qv (label long)
SCLIM = ["dycore_SC_restart_qv", "dycore_SC_output30_qv"]
SCLONG = [f"dycore_SC_restart_{f}" for f in "theta_m qv qc qr".split()] + \
    [f"dycore_SC_output30_{f}" for f in "qc qr".split()]
MUTATIONS += [
    # atm_advance_scalars_work (rk_steps 1, 2)
    (AS, "        sw2(k, j) = m.adv_coefs(j, iEdge) + sign1(uhAvg(k, iEdge)) * m.adv_coefs_3rd(j, iEdge);",
     "        sw2(k, j) = m.adv_coefs(j, iEdge) + m.adv_coefs_3rd(j, iEdge);", SCQ([1])),
    (AS, "sw2(k, 10) * scalar_new(iScalar, k, ica[10]);", "sw2(k, 10) * scalar_new(iScalar, k, ica[9]);", SCQ([1])),
    (AS, "          STC(iScalar, k) = STC(iScalar, k) - m.edgesOnCell_sign(i, iCell)",
     "          STC(iScalar, k) = STC(iScalar, k) + m.edgesOnCell_sign(i, iCell)", SCQ([1])),
    (AS, "      WDTN(iScalar, 2) = wwAvg(2, iCell) * (fnm(2) * scalar_new(iScalar, 2, iCell) + fnp(2) * scalar_new(iScalar, 2 - 1, iCell));",
     "      WDTN(iScalar, 2) = wwAvg(2, iCell) * (fnm(2) * scalar_new(iScalar, 2, iCell) + fnp(2) * scalar_new(iScalar, 2, iCell));",
     SCQ([1])),
    (AS, "                                 coef_3rd_order);\n    for (int k = 1; k <= nVertLevels; ++k)  // :4920-4928",
     "                                 1.0);\n    for (int k = 1; k <= nVertLevels; ++k)  // :4920-4928", SCQ([1])),
    (AS, "    weight_time_new = 1.;                     // :4750",
     "    weight_time_new = 0.5;                    // :4750", SCQ([1])),
    (AS, "    atm_advance_scalars_work(cfg, m, s, f, rk_timestep[rk_step], horiz_flux_array, rk_step,",
     "    atm_advance_scalars_work(cfg, m, s, f, rk_timestep[3], horiz_flux_array, rk_step,", SCQ([1])),
    # atm_advance_scalars_mono_work (rk_step 3)
    (AS, "          const int ii = uhAvg(k, iEdge) > 0 ? 1 : 2;", "          const int ii = uhAvg(k, iEdge) > 0 ? 2 : 1;", SCQ([1])),
    (AS, "        for (int j = 1; j <= 6; ++j) {", "        for (int j = 1; j <= 5; ++j) {", SCLIM),
    (AS, "        flux_upwind_arr[k] = dt * (fmax0(wwAvg(k, iCell)) * scalar_old(k - 1, iCell) +",
     "        flux_upwind_arr[k] = dt * (fmin0(wwAvg(k, iCell)) * scalar_old(k - 1, iCell) +", SCQ([1])),
    (AS, "        scale_arr(k, SCALE_IN, iCell) = -(rdnw(k) * (fmin0(wdtn(k + 1, iCell)) - fmax0(wdtn(k, iCell))));",
     "        scale_arr(k, SCALE_IN, iCell) = -(rdnw(k) * (fmin0(wdtn(k, iCell)) - fmax0(wdtn(k, iCell))));", SCQ([1])),
    # min(1, .) of SCALE_IN deleted is exact (shown on the 30-min state): every use of SCALE_IN is a min with a
    # SCALE_OUT that is <= 1 (TI:5680-5681, 5703-5704), and symmetrically for SCALE_OUT
    (AS, "        real scale_factor = (s_max(k, iCell) * rho_lim(k, iCell) - scalar_new(k, iCell)) /",
     "        real scale_factor = (s_max(k, iCell) * rho_zz_old(k, iCell) - scalar_new(k, iCell)) /", SCLIM),
    (AS, "          flux = fmax0(flux) * vmin(scale_arr(k, SCALE_OUT, cell1), scale_arr(k, SCALE_IN, cell2)) +",
     "          flux = fmax0(flux) * vmin(scale_arr(k, SCALE_OUT, cell1), scale_arr(k, SCALE_OUT, cell2)) +", SCQ([1, 2])),
    (AS, "        flux = fmax0(flux) * vmin(scale_arr(k - 1, SCALE_OUT, iCell), scale_arr(k, SCALE_IN, iCell)) +",
     "        flux = fmax0(flux) * vmin(scale_arr(k, SCALE_OUT, iCell), scale_arr(k, SCALE_IN, iCell)) +", SCQ([1, 2])),
    (AS, "          flux_arr(k, iEdge) = dt * flux_arr(k, iEdge) - flux_upwind;",
     "          flux_arr(k, iEdge) = dt * flux_arr(k, iEdge) + flux_upwind;", SCQ([1])),
    (AS, "      for (int k = 1; k <= nVertLevels; ++k) scalars_new(iScalar, k, iCell) = fmax0(scalar_new(k, iCell));",
     "      for (int k = 1; k <= nVertLevels; ++k) scalars_new(iScalar, k, iCell) = scalar_new(k, iCell);", SCLIM),
    (AS, "            (scalar_new(k, iCell) + (-(rdnw(k) * (wdtn(k + 1, iCell) - wdtn(k, iCell))))) / rho_lim(k, iCell);",
     "            (scalar_new(k, iCell) + (-(rdnw(k) * (wdtn(k + 1, iCell) - wdtn(k, iCell))))) / rho_zz_old(k, iCell);",
     SCQ([1])),
    # the call (TI:2504-2507) and the transport's checks
    (DY, "        advance_scalars(cfg, m, s, f, rk_step, rk_timestep);", "        (void)rk_timestep;", SCQ([1])),
    # the pentagon branches (TI:4801-4821, 5450-5463; RG's x1.10242 mesh has 12 pentagons)
    (AS, "            const real scalar_weight = m.adv_coefs(j, iEdge) + sign1(uhAvg(k, iEdge)) * m.adv_coefs_3rd(j, iEdge);",
     "            const real scalar_weight = m.adv_coefs(j, iEdge) - sign1(uhAvg(k, iEdge)) * m.adv_coefs_3rd(j, iEdge);",
     ["dycore_RG_step1_qv"]),
    (AS, "                  uhAvg(k, iEdge) * (m.adv_coefs(i, iEdge) + sign1(uhAvg(k, iEdge)) * m.adv_coefs_3rd(i, iEdge));",
     "                  uhAvg(k, iEdge) * (m.adv_coefs(i, iEdge) + m.adv_coefs_3rd(i, iEdge));", ["dycore_RG_step1_qv"]),
    (AS, "  if (!cfg.config_monotonic) no(", "  if (false) no(", ["unit_dycore"]),
    (DY, "  if (s.scalar_names.size() > 1 && !c.config_scalar_advection) no(", "  if (false) no(", ["unit_dycore"]),
    (DY, "  if (s.scalar_names != want)\n", "  if (false)\n", ["unit_dycore"]),
    # vertical mixing of the perturbation (DM:850-882, 1442-1466)
    (DP, "      u_mix[k] = u(k, iEdge) - u_init(k) * cosa + v_init(k) * sina;",
     "      u_mix[k] = u(k, iEdge) - u_init(k) * sina + v_init(k) * cosa;", SCT([1], "u")),
    (DP, "              (((theta_m(k + 1, iCell) - t_init(k + 1, iCell)) - (theta_m(k, iCell) - t_init(k, iCell))) / (zp - z0) -",
     "              (((theta_m(k + 1, iCell) - t_init(k, iCell)) - (theta_m(k, iCell) - t_init(k, iCell))) / (zp - z0) -",
     SCT([1], "theta_m")),
    (DS, "      u_vertical_mixing_perturbation(m, a, u, rho_edge, tend_u_euler);\n      return;",
     "      return;", SCT([10], "u") + ["dycore_SC_restart_u"]),
    (DS, "      theta_vertical_mixing_perturbation(m, a, theta_m, rho_zz, tend_theta_euler);\n      return;",
     "      return;", SCT([1], "theta_m")),
    # Kessler microphysics and its interface
    (DY, "    driver_microphysics(m, s, f.theta_m2, f.rho_zz2, f.scalars2, d, f.rt_diabatic_tend, f.mp);", "",
     SCT([1], "theta_m")),
    (MP, "      rt_diabatic_tend(k, i) = (theta_m(k, i) - rt_diabatic_tend(k, i)) / dt_dyn;",
     "      rt_diabatic_tend(k, i) = 0.0;", SCT([2], "theta_m")),
    (MP, "      exner(k, i) = libm::pow(zz(k, i) * (R_d / P0) * (rtheta_p(k, i) + rtheta_b(k, i)), rcv);",
     "      exner(k, i) = libm::mvec_pow(zz(k, i) * (R_d / P0) * (rtheta_p(k, i) + rtheta_b(k, i)), rcv);",
     SCT([2], "theta_m")),
    (MP, "      P(th_p, i, k) = theta_m(k, i) / (1.0 + R_v / R_d * fmax0(P(qv_p, i, k)));",
     "      P(th_p, i, k) = theta_m(k, i) / (1.0 + R_v / R_d * 0.5 * fmax0(P(qv_p, i, k)));", SCT([1], "theta_m")),
    (MP, "      theta_m(k, i) = P(th_p, i, k) * (1.0 + rvord * P(qv_p, i, k));",
     "      theta_m(k, i) = P(th_p, i, k) * (1.0 + rvord * 0.5 * P(qv_p, i, k));", SCT([1], "theta_m")),
    (MP, "      rtheta_p(k, i) = rho_zz(k, i) * theta_m(k, i) - rtheta_b(k, i);",
     "      rtheta_p(k, i) = rho_zz(k, i) * theta_m(k, i);", SCT([2], "theta_m")),
    (MP, "      pressure_p(k, i) = zz(k, i) * R_d * (exner(k, i) * rtheta_p(k, i) + (exner(k, i) - exner_b(k, i)) * rtheta_b(k, i));",
     "      pressure_p(k, i) = zz(k, i) * R_d * (exner(k, i) * rtheta_p(k, i) + (exner(k, i) - exner_b(k, i)) * rtheta_p(k, i));",
     SCT([1], "pressure_p")),
    (MP, "      scalars(index_qc, k, i) = P(qc_p, i, k);", "      scalars(index_qc, k, i) = P(qc_p, i, k) + 1.0e-12;",
     SCT(range(1, 11), "qc")),
    (MP, "      scalars(index_qr, k, i) = P(qr_p, i, k);", "      scalars(index_qr, k, i) = P(qr_p, i, k) + 1.0e-12;",
     SCT(range(1, 11), "qr")),
    (MP, "  mp.dt_microp = mp.dt_dyn;     // MGR:708", "  mp.dt_microp = 2.0 * mp.dt_dyn;  // MGR:708", SCLONG),
    (KS, "CPPMPAS_HD inline bool lane(int i, int n) { return n > 1 && !((n & 1) && i == n - 1); }",
     "CPPMPAS_HD inline bool lane(int i, int n) { return false && n > 1 && !((n & 1) && i == n - 1); }", SCLONG),
    (KS, "CPPMPAS_HD inline double rpow(double a, double b) { return ::cppmpas::libm::pow(a, b); }",
     "CPPMPAS_HD inline double rpow(double a, double b) { return ::cppmpas::libm::mvec_pow(a, b); }", SCLONG),
    (MP, "  if (mp.mp_top_level != nVertLevels)", "  if (false)", ["unit_microphysics"]),
    (MP, "    if (layer_height <= config_microphysics_top) mp.mp_top_level = k;", "", ["unit_microphysics"]),
    (KS, "      if (es_cap) es = min(es, real(0.99) * pressure);", "", ["unit_microphysics"]),
    # dropped as equivalent on SC: "kk.es_cap = true" -> false (es <= 0.99*pressure in every SC column, so the cap
    # never acts; unit_microphysics shows the cap itself); the source update of TI:5208 deleted (scalar_tend is 0,
    # it only turns a -0.0 of time level 1 into +0.0, and no -0.0 reaches it); fmax0(qv) -> qv in th_p (qv >= +0);
    # sincos -> sin / cos of angleEdge in u_mix (the SC angleEdge values pi/6, pi/2, 5pi/6 give the same bits); the
    # v_init term of u_mix deleted (SC's v_init is 0); R_v / R_d -> rvord in th_p (the same constant, 461.6/287);
    # theta_m from scalars(index_qv) instead of qv_p (equal: the scalar was set from qv_p one line above); the
    # first four were run and survive every SC step test, the last three survived the step-1 tests they had.
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
MUTATIONS += [
    # --- physics column checks, appended after the supercell's mutants (merge of the physics branches: the
    # supercell keeps its indices 493-534, main 535-538; WSM6 539-541, YSU / surface layer 542-551, radiation 552-564, Noah 565-573,
    # New Tiedtke 574-579, GWDO 580-582)
    # WSM6 column check (tests/wsm6_column): the libmvec lanes of MPAS's -O3 build (scalar log / exp everywhere,
    # CPPWRF's kernel as WRF calls it), the max(0, qv) of th_p (PI:665; cells of regime 4 have qv < 0), the
    # background radius above mp_top_level (PI:1010; case 9:41:2:33)
    ("src/physics/module_mp_wsm6.cpp", "pii.data(), p.data(), delz.data(), rain, rainncv, snow, snowncv, graupel, graupelncv, sr, i - its,\n           ite - its + 1);",
     "pii.data(), p.data(), delz.data(), rain, rainncv, snow, snowncv, graupel, graupelncv, sr, i - its,\n           1);", ["wsm6_column"]),
    ("src/physics/mpas_atmphys_microphysics.cpp", "P(th_p, i, k) = m.theta_m(k, i) / (1. + R_v / R_d * fmax0(P(qv_p, i, k)));",
     "P(th_p, i, k) = m.theta_m(k, i) / (1. + R_v / R_d * P(qv_p, i, k));", ["wsm6_column"]),
    ("src/physics/mpas_atmphys_microphysics.cpp", "        m.re_snow(k, i) = re_qs_bg;", "        m.re_snow(k, i) = re_qi_bg;", ["wsm6_column"]),
    # surface layer / YSU column check (tests/pbl_column): loop 30's thcon lanes (sf_sfclayrev.F90:258) off, WRF's
    # scalar table loop instead of MPAS's -O3 one (sf_sfclayrev_init), WRF 4.8.0's ent_eff (bl_ysu.F90:855; cloud-topped
    # columns with config_ysu_pblmix), MPAS_to_physics' top pressure without libmvec lanes (PI:504), qv_p without
    # max(0, qv) (PI:330; regime 5 has qv < 0 and -0.0), kzq from kzm (PB:604 takes kzh), the sea-ice regime kept
    # from xice >= 0.6 (SL:731; xice in [0.3, 1]), the sea-ice skin temperature floor (SL:524), the cell2 north term
    # of tend_toEdges (TD:547-549), the moist theta coupling (TD:492)
    ("src/physics/mpas_atmphys_sfclayer.cpp", "  k.thcon_lanes = true;", "  k.thcon_lanes = false;", ["pbl_column"]),
    ("src/physics/mpas_atmphys_sfclayer.cpp", "    sfclayrev::init_sfclayrev_tables_o3(tables);", "    sfclayrev::init_sfclayrev_tables(tables);", ["pbl_column"]),
    ("src/physics/mpas_atmphys_pbl.cpp", "    k.ent_eff_mpas = true;", "    k.ent_eff_mpas = false;", ["pbl_column"]),
    ("src/physics/mpas_atmphys_interface.cpp", "      if (lane(i - its, n))", "      if (false)", ["pbl_column"]),
    ("src/physics/mpas_atmphys_interface.cpp", "      v.qv_p(i, k) = fmax0(q(m.index_qv, k, i));", "      v.qv_p(i, k) = q(m.index_qv, k, i);", ["pbl_column"]),
    ("src/physics/mpas_atmphys_pbl.cpp", "      m.kzq(k, i) = kzh_p(i, k);", "      m.kzq(k, i) = kzm_p(i, k);", ["pbl_column"]),
    ("src/physics/mpas_atmphys_sfclayer.cpp", "        if (x >= 0.5) m.regime(i)", "        if (x >= 0.6) m.regime(i)", ["pbl_column"]),
    ("src/physics/mpas_atmphys_sfclayer.cpp", "        sea.tsk[I(i)] = sfclayrev::c_max(sst_p[I(i)], real(271.4));", "        sea.tsk[I(i)] = sst_p[I(i)];", ["pbl_column"]),
    ("src/physics/mpas_atmphys_todynamics.cpp", "+ Uy_tend(k, cell2) * 0.5 * (n(1, iEdge) * m.north(1, cell2)", "+ Uy_tend(k, cell2) * 0.5 * (n(1, iEdge) * m.east(1, cell2)", ["pbl_column"]),
    ("src/physics/mpas_atmphys_todynamics.cpp", "      tend_th(k, i) = coeff * tend_th(k, i) + R_v / R_d", "      tend_th(k, i) = tend_th(k, i) + R_v / R_d", ["pbl_column"]),
    # --- radiation (tests/radiation_column): cloudiness, RRTMG wrappers and drivers, ozone climatology
    # LW buffer layer count (LW 11822)
    (RAD, 'std::max(static_cast<int>(nint(pw1d_top / deltap)), 1)',
     'std::max(static_cast<int>(nint(pw1d_top / deltap)), 2)', ["radiation_column"]),
    # table temperature profile matched at level kte (LW 11993)
    (RAD, 'v.tz[K(k)] = w.varint[K(k)] + (v.tz[K(kte)] - w.varint[K(kte)]);',
     'v.tz[K(k)] = w.varint[K(k)] + (v.tz[K(kte + 1)] - w.varint[K(kte + 1)]);', ["radiation_column"]),
    # vinterp_ozn interpolation weights (VI:126-127)
    (RADH, '        V(k, (O(kupper) * dpl + O(kupper + 1) * dpu) / (dpl + dpu));\n        found = true;',
     '        V(k, (O(kupper) * dpu + O(kupper + 1) * dpl) / (dpl + dpu));\n        found = true;', ["radiation_column"]),
    # LW inirad ozone on every layer with config_o3climatology = false (LW 12024-12026)
    (RAD, 'for (int k = kts; k <= nlayers; ++k) v.wkl[2][K(k)] = v.o3mmr[K(k)] * amdo;',
     'for (int k = kts; k <= nlayers - 1; ++k) v.wkl[2][K(k)] = v.o3mmr[K(k)] * amdo;', ["radiation_column"]),
    # relcalc snow ramp (LW 12493)
    (RAD, '    r = r + (rliqocean - r) * std::min(1.0, std::max(0.0, snowh * 10.));',
     '    r = r + (rliqocean - r) * std::min(1.0, std::max(0.0, snowh * 1.));', ["radiation_column"]),
    # LW McICA permuteseed 150 (LW 12183)
    (RAD, 'mcica_subcol(icld, 150, 140,',
     'mcica_subcol(icld, 151, 140,', ["radiation_column"]),
    # zenith angle: the reference's sincos of xxlat (SW 10205). SURVIVES: on the check's latitudes sincos and sin / cos
    # round the same (s_sincos.c differs from __sin only in some ulps of 0.855 <= |x| < 2.43); kept as a known survivor
    (RAD, '  libm::sincos(xxlat, &sinlat, &coslat);',
     '  sinlat = libm::sin(xxlat); coslat = libm::cos(xxlat);', ["radiation_column"]),
    # SW flux adjustment scon / rrsw_scon (inatm_sw)
    (RAD, '  const rb rrsw_scon = rb(1.36822e+03);',
     '  const rb rrsw_scon = rb(1.3682e+03);', ["radiation_column"]),
    # radconst eccentricity factor (RS:1045)
    (RDRV, '  const real eccfac = 1.000110 + 0.034221 * c1 + 0.001280 * s1 +',
     '  const real eccfac = 1.000110 + 0.034221 * c1 + 0.001280 * c1 +', ["radiation_column"]),
    # calc_cldfraction exp argument floor (CL:358)
    (RDRV, '        const real arg = std::max(-6.9, -alpha0 * qcld / denom);',
     '        const real arg = std::max(-6.8, -alpha0 * qcld / denom);', ["radiation_column"]),
    # o3climatology_from_MPAS day offset (OZ:279)
    (RDRV, '  double intjulian_r8 = double(julian + 1.0);  // OZ:279',
     '  double intjulian_r8 = double(julian);  // OZ:279', ["radiation_column"]),
    # init_o3climatology latitude interpolation order (OZ:178-179)
    (RDRV, '            ozmixm(mm, k, iCell) = OZI(j, k, i1, mm) + (OZI(j, k, i2, mm) - OZI(j, k, i1, mm)) * dlatCell / dlat;',
     '            ozmixm(mm, k, iCell) = OZI(j, k, i1, mm) + (OZI(j, k, i2, mm) - OZI(j, k, i1, mm)) * (dlatCell / dlat);', ["radiation_column"]),
    # o3vmr diagnostic pressures in hPa (RL:489)
    (RDRV, '        p2dv(i, k) = pres_hyd_p(i, k) / 100.;',
     '        p2dv(i, k) = pres_hyd_p(i, k) / 10.;', ["radiation_column"]),
    # Noah column check (tests/noah_column): log10 of TDFCND, the soil table read (SOILPARM.TBL SATDK / SATDW),
    # lsminit's supercooled water with MPAS's gravity (WRF's lsminit uses 9.81), the convective rain of
    # lsm_from_MPAS (:460-468), correct_tsk_over_seaice's xice < 0.1 limit (lsm_shared.F:49), usemonalb =
    # config_sfc_albedo (driver_lsm.F:735, case 'as'), the sfc-diagnostics cqs2 threshold (regime 9), the sea-ice
    # TBOT (seaice_drv.F, shared kernel), SFLX_GLACIAL's LVCOEF (snow ageing on land ice)
    (L, "CPPMPAS_FORCEINLINE double log10(double x) { return ::cppwrf::glibc::log10(arg(x)); }",
     "CPPMPAS_FORCEINLINE double log10(double x) { return ::cppwrf::glibc::log(arg(x)) / 2.302585092994046; }",
     ["noah_column"]),
    ("src/physics/mpas_atmphys_lsm_noahinit.cpp", "&t.satpsi[lc], &t.satdk[lc], &t.satdw[lc], &t.wltsmc[lc]",
     "&t.satpsi[lc], &t.satdw[lc], &t.satdk[lc], &t.wltsmc[lc]", ["noah_column"]),
    ("src/physics/mpas_atmphys_driver_lsm.cpp", "real fk = libm::pow((hlice / (pc::gravity * (-psisat)))",
     "real fk = libm::pow((hlice / (9.81 * (-psisat)))", ["noah_column"]),
    ("src/physics/mpas_atmphys_driver_lsm.cpp", '  if (cfg.config_convection_scheme != "off")  // :460-468',
     "  if (false)  // :460-468", ["noah_column"]),
    ("src/physics/mpas_atmphys_driver_lsm.cpp", "      if (xice[ii] < 0.1 && tsk_ice[ii] < 263.15) tsk_ice[ii] = 263.15;",
     "      if (false && tsk_ice[ii] < 263.15) tsk_ice[ii] = 263.15;", ["noah_column"]),
    ("src/physics/mpas_atmphys_driver_lsm.cpp", "a.usemonalb = cfg.config_sfc_albedo;", "a.usemonalb = true;",
     ["noah_column"]),
    ("src/physics/mpas_atmphys_driver_lsm.cpp", "    if (p.cqs2(i) < 1.e-5) {", "    if (p.cqs2(i) < 1.e-6) {",
     ["noah_column"]),
    ("src/physics/module_sf_noahdrv_kernels.hpp", "    real tbot = real(271.36);", "    real tbot = real(271.15);",
     ["noah_column"]),
    ("src/physics/module_sf_noahdrv.cpp", "gl.lvcoef = tables.lvcoef_data;", "gl.lvcoef = tables.salp_data;",
     ["noah_column"]),
    # New Tiedtke (cu_ntiedtke) and its MPAS driver, column check (tests/cu_gwdo_column): 574 cuinin's cuadjtqn
    # without the -O3 libmvec exp lanes (CPPWRF / WRF's scalar kernel), 575 dx_p (CV:840) without the
    # _ZGVbN2vv_pow lanes, 576 initflag without config_do_restart (CV:478; restart cases), 577 thften without the
    # shortwave heating (CV:846-847), 578 pratec = rn / dt without stepcu (PO:96; the n_cu = 2 case), 579 the
    # rainc bucket without the i_rainc count (CV:1094; acrain cases)
    ("src/physics/module_cu_ntiedtke.cpp", "    c.cuinin_mvec = lane(i - its, im);", "    c.cuinin_mvec = false;", ["cu_gwdo_column"]),
    ("src/physics/mpas_atmphys_convection.cpp", "    dx_p(i) = m.config_len_disp / (lane(i - its, ite - its + 1) ? libm::mvec_pow", "    dx_p(i) = m.config_len_disp / (false ? libm::mvec_pow", ["cu_gwdo_column"]),
    ("src/physics/mpas_atmphys_convection.cpp", "  if (m.config_do_restart || itimestep > 1) initflag = 0;", "  if (itimestep > 1) initflag = 0;", ["cu_gwdo_column"]),
    ("src/physics/mpas_atmphys_convection.cpp", "          (m.rthdynten(k, i) + m.rthblten(k, i) + m.rthratenlw(k, i) + m.rthratensw(k, i)) * v.pi_p(i, k);", "          (m.rthdynten(k, i) + m.rthblten(k, i) + m.rthratenlw(k, i)) * v.pi_p(i, k);", ["cu_gwdo_column"]),
    ("src/physics/module_cu_ntiedtke.cpp", "    (*a.pratec)(i) = rn / (real(a.stepcu) * a.dt);", "    (*a.pratec)(i) = rn / a.dt;", ["cu_gwdo_column"]),
    ("src/physics/mpas_atmphys_convection.cpp", "      m.i_rainc(i) = m.i_rainc(i) + 1;", "      m.i_rainc(i) = m.i_rainc(i);", ["cu_gwdo_column"]),
    # GWDO (bl_ysu_gwdo) and its MPAS driver, column check (tests/cu_gwdo_column): 580 WRF 4.8.0's cg = 1 (MPAS's
    # MMM-physics a4baf7f bl_gwdo.F90:153 has 0.5), 581 WRF 4.8.0's olmin floor on ol (GW:342 has none; regime 4
    # has ol = 0), 582 rubldiff without the PBL tendency removed (GD:680)
    ("src/physics/bl_gwdo.hpp", "const T xl = T(4.0e4), gmax = T(1.), veleps = T(1.0), frc = T(1.0), ce = T(0.8), cg = T(0.5);", "const T xl = T(4.0e4), gmax = T(1.), veleps = T(1.0), frc = T(1.0), ce = T(0.8), cg = T(1.);", ["cu_gwdo_column"]),
    ("src/physics/bl_gwdo.hpp", "    ol = ol4[(nwd - 1) % 4 + 1];", "    ol = c_max(ol4[(nwd - 1) % 4 + 1], olmin);", ["cu_gwdo_column"]),
    ("src/physics/mpas_atmphys_gwdo.cpp", "      m.rubldiff(k, i) = rublten_p(i, k) - m.rublten(k, i);", "      m.rubldiff(k, i) = rublten_p(i, k);", ["cu_gwdo_column"]),
]

# ---------------------------------------------------------------------------------------------------------------
# the mesoscale_reference suite in the model run against RG (mpas842-dbl-rg-004): the split dynamics-transport, the
# physics manager / init / driver and their coupling (tests dycore_RG_*). Indices 583-599.
DRV, LU = "src/physics/mpas_atmphys_driver.cpp", "src/physics/mpas_atmphys_landuse.cpp"
RG1 = lambda *fs: [f"dycore_RG_step1_{f}" for f in fs]
MUTATIONS += [
    # split dynamics-transport: dt_dynamics (TI:2055), the transport's rk_timestep(1) of order 2 (TI:2679), the
    # weights of the density (TI:4754), rho_zz_int (TI:5270), the ruAvg sum of the substeps (TI:7790), the time level
    # copy between substeps (TI:7752)
    (DY, "    dt_dynamics = dt / static_cast<real>(dynamics_split);", "    dt_dynamics = dt;", RG1("u", "theta_m")),
    (DY, "    if (cfg.config_time_integration_order == 2) rk_timestep[1] = dt / 2.;", "", RG1("qv")),
    (AS, "    if (rk_step == 2) weight_time_new = 1. / 2;", "    if (rk_step == 2) weight_time_new = 1.;", RG1("qv")),
    (AS, "            rho_zz_old(k, iCell) + dt * (rho_zz_int(k, iCell) - rdnw(k) * (wwAvg(k + 1, iCell) - wwAvg(k, iCell)));",
     "            rho_zz_old(k, iCell) + dt * rho_zz_int(k, iCell);", RG1("qv")),
    (DY, "f.ruAvg_split(k, iEdge) = f.ruAvg(k, iEdge) + f.ruAvg_split(k, iEdge);",
     "f.ruAvg_split(k, iEdge) = f.ruAvg(k, iEdge);", RG1("qv")),
    (DY, "        s.theta_m(k, iCell) = f.theta_m2(k, iCell);\n", "\n", RG1("theta_m")),
    # physics_run_init / physics_init / landuse_init_forMPAS / physics_update_surface / physics_timetracker
    (DRV, "  r.xice_threshold = c.config_frac_seaice ? 0.02 : 0.5;", "  r.xice_threshold = c.config_frac_seaice ? 0.5 : 0.02;",
     RG1("theta_m")),
    (DRV, "      dzs(1, iCell) = 0.10;", "      dzs(1, iCell) = 0.11;", RG1("theta_m")),
    (DRV, "      east(1, iCell) = -libm::sin(lon);", "      east(1, iCell) = libm::sin(lon);", RG1("u")),
    # SURVIVES (shown once): the greenness update at step 1 recomputes the init file's own sfc_albbck / vegfra / shdmin
    # / shdmax for the same date (init_atmosphere ran the same monthly_interp_to_date), and the next ring is at 24 h,
    # after the last output of every reference; the update itself is killed by the shdmin mutant below
    (DRV, "  if (r.green_alarm.ringing(now_s)) {", "  if (false && r.green_alarm.ringing(now_s)) {", RG1("theta_m")),
    (DRV, "      o3climatology_from_MPAS(r.rad.curr_julday, q, m.nCellsSolve());",
     "      o3climatology_from_MPAS(r.rad.curr_julday + 1.0, q, m.nCellsSolve());", RG1("theta_m")),
    ("src/physics/mpas_atmphys_driver.hpp", "    while (ring <= now) ring += interval;", "    while (ring < now) ring += interval;",
     ["dycore_RG_step2_theta_m"]),
    (DRV, "      rqvdynten(k, i) = config_monotonic ? (scalars2(qv, k, i) - s.scalars(qv, k, i)) / dt : 0.;",
     "      rqvdynten(k, i) = 0.;", ["dycore_RG_step2_theta_m"]),
    (LU, "    if (julday < 105 || julday >= 288) isn = 2;", "    if (julday < 105 || julday >= 285) isn = 2;", RG1("theta_m")),
    # SURVIVES (shown once): over water sfc_albbck is read by no scheme after the update (the open-water albedo comes
    # from landuse_init_forMPAS's sfc_albedo at the start), so the 0.08 is exact on RG
    (LU, "      if (landmask(iCell) == 0) sfc_albbck(iCell) = 0.08;", "      if (landmask(iCell) == 0) sfc_albbck(iCell) = 0.07;",
     RG1("theta_m")),
    (LU, "        field_out(n) = (field_in(month2, n) * real(target_date - middle[l]) +",
     "        field_out(n) = (field_in(month1, n) * real(target_date - middle[l]) +", RG1("theta_m")),
    # monthly_min_max's minimum (shdmin, read by Noah) in the step-1 greenness update (DT:228-258)
    (LU, "      field_min(n) = minner;", "      field_min(n) = maxxer;", RG1("theta_m")),
    # the writer keeps the input file's attributes of the physics fields (the reference restart's isltyp)
    (W, "      if (std::find(rd.begin(), rd.end(), n) != rd.end()) {", "      if (false) {", ["dycore_RG_restart_isltyp"]),
    # WSM6 output diagnostics (l_diags): the wsm6_column cases with l_diags (refl10cm*, relhum)
    (RF, "      ilamr[k] = 1. / lamr;", "      ilamr[k] = 1.01 / lamr;", ["wsm6_column"]),
    (RF, "        const double fmelt_s = c_max(0.005, c_min(1.0 - rs[k] / rs[k_0], 0.99));",
     "        const double fmelt_s = c_max(0.05, c_min(1.0 - rs[k] / rs[k_0], 0.99));", ["wsm6_column"]),
    (RF, "        const double fmelt_g = c_max(0.005, c_min(1.0 - rg[k] / rg[k_0], 0.99));",
     "        const double fmelt_g = c_max(0.05, c_min(1.0 - rg[k] / rg[k_0], 0.99));", ["wsm6_column"]),
    (RF, "    dBZ[k] = 10. * (lane ? libm::mvec_log10(ze) : libm::log10(ze));", "    dBZ[k] = 10. * libm::log10(ze);",
     ["wsm6_column"]),
    (RD, "  tmp = (x + 0.5) * ::cppmpas::libm::mvec_log(tmp) - tmp;", "  tmp = (x + 0.5) * ::cppmpas::libm::log(tmp) - tmp;",
     ["wsm6_column"]),
    (RD, "  mra = mra + (1.0 - mra) * fm;", "  mra = mra + (1.0 - mra) * fm * 0.5;", ["wsm6_column"]),
    (MC, "        dBZ1d[k - 1] = -35. > dBZ1d[k - 1] ? -35. : dBZ1d[k - 1];  // max(-35._RKIND, dBZ1d(k))",
     "        dBZ1d[k - 1] = -30. > dBZ1d[k - 1] ? -30. : dBZ1d[k - 1];  // max(-35._RKIND, dBZ1d(k))", ["wsm6_column"]),
    (MC, "        if (zp[k - 1] < 1000.) kp = k;", "        if (zp[k - 1] < 1100.) kp = k;", ["wsm6_column"]),
    (MC, "        if (dBZ1d[k - 1] > mx) mx = dBZ1d[k - 1];", "        if (dBZ1d[k - 1] < mx) mx = dBZ1d[k - 1];",
     ["wsm6_column"]),
    (MC, "      m.refl10cm_1km_max(i) = m.refl10cm_1km_max(i) > m.refl10cm_1km(i) ? m.refl10cm_1km_max(i) : m.refl10cm_1km(i);",
     "      m.refl10cm_1km_max(i) = m.refl10cm_1km(i);", ["wsm6_column"]),
    (MC, "        if (tempc <= 0.) qvs1 = rsif(p1, t1);", "        if (tempc <= -1.) qvs1 = rsif(p1, t1);", ["wsm6_column"]),
    (MC, "                 C8 = -.321582393E-13;", "                 C8 = -.321582393E-12;", ["wsm6_column"]),
    (MC, "                 C8 = .161444444E-12;", "                 C8 = .161444444E-11;", ["wsm6_column"]),
    (MC, "  const real X = xt > -80. ? xt : -80.;\n  real ESI", "  const real X = xt > -60. ? xt : -60.;\n  real ESI",
     ["wsm6_column"]),
    (MC, "        m.relhum(k, i) = m.relhum(k, i) * 100.;", "        m.relhum(k, i) = m.relhum(k, i) * 99.;", ["wsm6_column"]),
]

# isobaric diagnostics (core_atmosphere/diagnostics/mpas_isobaric_diagnostics.cpp, mpas_isobaric_diagnostics.F)
ISO = "src/core_atmosphere/diagnostics/mpas_isobaric_diagnostics.cpp"
I0 = lambda *fs: [f"diag_RG_iso0_{f}" for f in fs]
MUTATIONS += [
    # interp_tofixed_pressure (IS:1027-1110): the weights of the bracketing levels, the below-ground value
    (ISO, "          const real dpu = pres_out(icol, k) - pres_in(icol, kupper[icol]);",
     "          const real dpu = pres_out(icol, k) - pres_in(icol, kupper[icol] + 1);", I0("temperature_500hPa")),
    (ISO, "        field_out(icol, k) = field_in(icol, nlev_in);", "        field_out(icol, k) = field_in(icol, nlev_in - 1);",
     I0("temperature_925hPa")),
    # pressure2 (IS:654-688): the middle levels, the lowest level, the libmvec lanes of the three loops
    (ISO, "      const real w2 = (height(k + 1, iCell) - height(k, iCell)) / (height(k + 1, iCell) - height(k - 1, iCell));",
     "      const real w2 = (height(k + 1, iCell) - height(k, iCell)) / (height(k + 1, iCell) - height(k, iCell));",
     I0("height_500hPa")),
    (ISO, "      const real z2 = 0.5 * (height(k + 1, iCell) + height(k + 2, iCell));",
     "      const real z2 = 0.5 * (height(k + 1, iCell) + height(k + 1, iCell));", I0("height_925hPa")),
    (ISO, "  const int nPair = nCells - nCells % 2;", "  const int nPair = 0;", I0("height_500hPa")),
    # pressure at vertices (IS:691-700)
    (ISO, "      pressure_v(k, iVert) = pressure_v(k, iVert) / areaTriangle(iVert);",
     "      pressure_v(k, iVert) = pressure_v(k, iVert) / areaTriangle(1);", I0("vorticity_500hPa")),
    # temperature and dewpoint (IS:702-716), the libmvec lanes of the dewpoint loop
    (ISO, "            (theta_m(k, iCell) / (1. + rvord * scalars(index_qv, k, iCell))) * exner(k, iCell);",
     "            (theta_m(k, iCell) / (1. + rvord * scalars(index_qv, k, iCell)));",
     I0("temperature_500hPa") + ["diag_RG24_iso3h_temperature_500hPa"]),
    (ISO, "        dewpoint(k, iCell) = dewpoint(k, iCell) + 273.15;", "        dewpoint(k, iCell) = dewpoint(k, iCell) + 273.16;",
     I0("dewpoint_500hPa")),
    (ISO, "    const int kPair = nVertLevels - nVertLevels % 2;", "    const int kPair = 0;", I0("dewpoint_850hPa")),
    # the interpolated sources (IS:767-902)
    (ISO, "    mass_level(need.relhum, relhum, \"relhum\");", "    mass_level(need.relhum, temperature, \"relhum\");",
     I0("relhum_500hPa")),
    (ISO, "    mass_level(need.uzonal, uzonal, \"uzonal\");", "    mass_level(need.uzonal, umeridional, \"uzonal\");",
     I0("uzonal_500hPa")),
    (ISO, "    w_level(vvel, \"w\");", "    w_level(height, \"w\");", I0("w_500hPa")),
    # need flags (IS:104-302): a group not computed
    (ISO, "  group(\"umeridional\", need.umeridional);", "  group(\"umeridional\", need.uzonal);", I0("umeridional_500hPa")),
    # compute_slp (IS:1113-1243): virtual temperature, the mm5 correction, t_surf's power, the libmvec lanes
    (ISO, "    const real thi = t(khi, icol) * (1. + 0.608 * scalars(index_qv, khi, icol));",
     "    const real thi = t(khi, icol) * (1. + 0.61 * scalars(index_qv, khi, icol));", I0("mslp")),
    (ISO, "        t_msl[icol] = tc - 0.005 * (d * d);", "        t_msl[icol] = tc - 0.005 * d;", I0("mslp")),
    (ISO, "    t_surf[icol] = t_at_pconst * libm::pow(p(1, icol) / p_at_pconst, gamma_rr_grav);  // pow@plt (0x23edc4)",
     "    t_surf[icol] = t_at_pconst * libm::pow(p(1, icol) / p_at_pconst, gamma * rr);", I0("mslp")),
    (ISO, "    const bool lane = icol <= ncol - ncol % 2;", "    const bool lane = false;", I0("mslp")),
    # compute_layer_mean (IS:1262-1351) and the isobaric levels (IS:616-638, 973-1046)
    (ISO, "        layerMean(iCell) = layerMean(iCell) + w * 0.5 * (field_in(iCell, k) + field_in(iCell, k + 1));",
     "        layerMean(iCell) = layerMean(iCell) + w * 0.5 * (field_in(iCell, k) + field_in(iCell, k));",
     I0("meanT_500_300")),
    (ISO, "    t_iso_levels(2) = 35000.0;", "    t_iso_levels(2) = 36000.0;", I0("t_iso_levels", "t_isobaric")),
    (ISO, "    for (int k = 1; k <= 13; ++k) z_iso_levels(k) = 30000.0 + 5000.0 * (k - 1);  // 30000.0 ... 90000.0, exact",
     "    for (int k = 1; k <= 13; ++k) z_iso_levels(k) = 30000.0 + 5000.0 * k;", I0("z_iso_levels", "z_isobaric")),
    (ISO, "        press_in(iCell, nVertLevelsP1 + 1 - k) = pressure2(k, iCell) * 100.0;",
     "        press_in(iCell, nVertLevelsP1 + 1 - k) = pressure2(k, iCell) * 100.1;", I0("z_isobaric")),
    # convective diagnostics (mpas_convective_diagnostics.F): t = 0 of RG (diag_RG_conv0_*) and 3 h / 6 h of RG24
    # (diag_RG24_conv_*, label long)
    (CV, "      dewpoint(k, iCell) = dewpoint(k, iCell) + 273.15;", "      dewpoint(k, iCell) = dewpoint(k, iCell) + 273.16;",
     CV0("dewpoint_surface")),
    (CV, "          (theta_m(k, iCell) / (1. + rvord * scalars(index_qv, k, iCell))) * exner(k, iCell);",
     "          (theta_m(k, iCell) / (1. + rvord * scalars(index_qv, k, iCell)));", CV0("temperature_surface")),
    (CV, "      const real lg = k <= nvec ? libm::mvec_log(evp / 6.112) : libm::log(evp / 6.112);",
     "      const real lg = false ? libm::mvec_log(evp / 6.112) : libm::log(evp / 6.112);", CV0("dewpoint_surface")),
    (CV, "    uzonal_surface(iCell) = uzonal(1, iCell);", "    uzonal_surface(iCell) = uzonal(2, iCell);", CV0("uzonal_surface")),
    (CV, "    umeridional_surface(iCell) = umeridional(1, iCell);", "    umeridional_surface(iCell) = umeridional(2, iCell);",
     CV0("umeridional_surface")),
    (CV, "    if (need_uzonal_1km) uzonal_1km(iCell) = column_height_value(ucol.data(), zp.data(), 1000.0, nVertLevels);",
     "    if (need_uzonal_1km) uzonal_1km(iCell) = column_height_value(ucol.data(), zp.data(), 1001.0, nVertLevels);",
     CV0("uzonal_1km")),
    (CV, "    if (need_umerid_1km) umeridional_1km(iCell) = column_height_value(vcol.data(), zp.data(), 1000.0, nVertLevels);",
     "    if (need_umerid_1km) umeridional_1km(iCell) = column_height_value(vcol.data(), zp.data(), 1001.0, nVertLevels);",
     CV0("umeridional_1km")),
    (CV, "    if (need_uzonal_6km) uzonal_6km(iCell) = column_height_value(ucol.data(), zp.data(), 6000.0, nVertLevels);",
     "    if (need_uzonal_6km) uzonal_6km(iCell) = column_height_value(ucol.data(), zp.data(), 6001.0, nVertLevels);",
     CV0("uzonal_6km")),
    (CV, "    if (need_umerid_6km) umeridional_6km(iCell) = column_height_value(vcol.data(), zp.data(), 6000.0, nVertLevels);",
     "    if (need_umerid_6km) umeridional_6km(iCell) = column_height_value(vcol.data(), zp.data(), 6001.0, nVertLevels);",
     CV0("umeridional_6km")),
    (CV, "  const real wzp1 = 1. - wz;", "  const real wzp1 = 1. + wz;", CV0("uzonal_1km umeridional_6km")),
    (CV, "      const real u_storm = u_mean + dev_motion * v_shear / shear_magnitude;",
     "      const real u_storm = u_mean - dev_motion * v_shear / shear_magnitude;", CV0("srh_0_1km srh_0_3km")),
    (CV, "      if (need_srh_01km) srh_0_1km(iCell) = integral_zpoint(srh.data(), zrel.data(), 0.0, 1000.0, nVertLevelsP1);",
     "      if (need_srh_01km) srh_0_1km(iCell) = integral_zpoint(srh.data(), zrel.data(), 0.0, 900.0, nVertLevelsP1);",
     CV0("srh_0_1km")),
    (CV, "      if (need_srh_03km) srh_0_3km(iCell) = integral_zpoint(srh.data(), zrel.data(), 0.0, 3000.0, nVertLevelsP1);",
     "      if (need_srh_03km) srh_0_3km(iCell) = integral_zpoint(srh.data(), zrel.data(), 0.0, 2900.0, nVertLevelsP1);",
     CV0("srh_0_3km")),
    (CV, "    const real midpoint_value = column_values[k - 1] + (column_values[k] - column_values[k - 1]) * zr_midpoint;",
     "    const real midpoint_value = column_values[k - 1];", CV0("srh_0_1km")),
    # (cape = cape + parea survives: every branch gives parea >= 0, so MAX(0.0, parea) is exact)
    (CV, "      parea = 0.5 * dz * (b1 + b2);", "      parea = 0.6 * dz * (b1 + b2);", CV0("cape")),
    (CV, "      cin = cin + narea;", "      cin = cin - narea;", CV0("cin")),
    (CV, "    const bool lane = k <= nvec;", "    const bool lane = false;", CV0("cape")),
    # getthe only picks the source parcel (kmax): 0.00125 -> 0.00126 in tlcl survives (no column changes its kmax)
    (CV, "         libm::exp(((3376.0 / tlcl) - 2.54) * q * (1.0 + 0.81 * q));",
     "         libm::exp(((3076.0 / tlcl) - 2.54) * q * (1.0 + 0.81 * q));", CV0("cape")),
    (CV, "        qv2 = minsd(qt, fliq * getqvs(p2, t2) + fice * getqvi(p2, t2));", "        qv2 = minsd(qt, 1.01 * getqvs(p2, t2));",
     CV0("cape")),
    (CV, "  (void)p.R(\"lcl\");", "  p.R(\"lcl\")(1) = 1.;", CV0("lcl")),
    (CV, "  (void)p.R(\"lfc\");", "  p.R(\"lfc\")(1) = 1.;", CV0("lfc")),
    (CV, "      const real uph = integral_zstaggered(&updraft_helicity(1, iCell), z_agl.data(), 2000.0, 5000.0, nVertLevels);",
     "      const real uph = integral_zstaggered(&updraft_helicity(1, iCell), z_agl.data(), 1000.0, 5000.0, nVertLevels);",
     CV24("updraft_helicity_max")),
    (CV, "          updraft_helicity(k, iCell) = updraft_helicity(k, iCell) + m.kiteAreasOnVertex(j, i) * vorticity(k, i);",
     "          updraft_helicity(k, iCell) = updraft_helicity(k, iCell) + m.kiteAreasOnVertex(j, i) * vorticity(k, i) * 1.5;",
     CV24("updraft_helicity_max")),
    (CV, "      w_velocity_max(iCell) = maxsd(mx, w_velocity_max(iCell));",
     "      w_velocity_max(iCell) = maxsd(0.5 * mx, w_velocity_max(iCell));", CV24("w_velocity_max")),
    (CV, "          maxsd(std::sqrt(uzonal(1, iCell) * uzonal(1, iCell) + umeridional(1, iCell) * umeridional(1, iCell)),",
     "          maxsd(std::sqrt(uzonal(1, iCell) * uzonal(2, iCell) + umeridional(1, iCell) * umeridional(1, iCell)),",
     CV0("wind_speed_level1_max")),
    (CV, "    if (c.field_will_be_written(f)) {\n      auto& a = c.pools->R(f);", "    if (false) {\n      auto& a = c.pools->R(f);",
     CV24("w_velocity_max", "06")),
]

# the output streams of cppmpas run --streams (io/mpas_streams.cpp, tools/cppmpas.cpp, the namelist attributes of
# io/mpas_writer.cpp): RG's forecast directory (dycore_RG_streams_*, about 50 s) and RG24's (dycore_RG24_*, long)
ST, NT = "src/io/mpas_streams.cpp", "src/core/namelist_table.inc"
RGS = lambda *fs: [f"dycore_RG_streams_{f}" for f in fs]
MUTATIONS += [
    # the output interval in seconds: every alarm shifted, the 01:00 history is not written
    (ST, "  return ((days * 24 + h) * 60 + m) * 60 + sec;", "  return ((days * 24 + h) * 60 + m) * 60 + sec + 600;",
     RGS("history.2026-10-12_01.00.00", "all")),
    # no restart file at the initial time (mpas_atm_core.F:667)
    (ST, '  if (OutputStream* r = sm.find("restart")) r->reset(0);', '  if (OutputStream* r = sm.find("restart_")) r->reset(0);',
     RGS("all")),
    # restart_timestamp as Fortran's list-directed write writes it
    (X, '            std::fprintf(tf, " %s\\n", stamp.c_str());', '            std::fprintf(tf, "%s\\n", stamp.c_str());', RGS("all")),
    # MG:615-645: the diagnostics alarm takes the shorter of the output and diagnostics intervals (relhum at 01:00)
    (X, "          step = std::min(step, interval_seconds(dia->spec.output_interval));",
     "          step = std::max(step, interval_seconds(dia->spec.output_interval));", RGS("history.2026-10-12_01.00.00")),
    # atm_reset_diagnostics (mpas_atm_core.F:884-895): refl10cm_1km_max restarts after each diagnostics write
    (X, "            for (std::size_t k = 0; k < r.size(); ++k) r.data()[k] = 0.;",
     "            for (std::size_t k = 0; k < r.size(); ++k) r.data()[k] = r.data()[k];",
     ["dycore_RG24_diag.2026-10-12_06.00.00"]),
    # the namelist options as global attributes: logicals, the schemes physics_namelist_check resolved
    (W, '        f.put_att_text("", name, yes_no(value == "true"));', '        f.put_att_text("", name, yes_no(value != "true"));',
     RGS("history.2026-10-12_00.00.00")),
    (X, '            {"config_microp_scheme", &pr.cfg.config_microp_scheme},',
     '            {"config_microp_scheme_", &pr.cfg.config_microp_scheme},', RGS("history.2026-10-12_00.00.00")),
    # a field read from the input stream keeps the input file's attribute values (dvEdge's units '-')
    (W, "            if (fa.name == k && fa.type == NC_CHAR) v = trimmed(fa.text());",
     "            if (fa.name == k && fa.type == NC_INT) v = trimmed(fa.text());", RGS("history.2026-10-12_00.00.00")),
    # an inactive field of a stream is skipped
    (ST, "        if (!packages.active(*rv)) continue;  // an inactive field is skipped",
     "        if (false) continue;  // an inactive field is skipped", ["dycore_run_RG_streams"]),
    # MPAS_field_will_be_written answered by the stream manager: the diagnostics stream's fields are computed
    (X, "      if (streams_mode) return sm.field_will_be_written(n, clock_s);", "      if (streams_mode) return false;",
     RGS("diag.2026-10-12_00.00.00")),
    # qv_init: the restart stream writes its Registry default
    (R, '    const bool mesh_default = rv && rv->name == "qv_init";', '    const bool mesh_default = rv && rv->name == "qv_init_";',
     ["dycore_run_RG_streams"]),
    # the RG24 forecast's config_run_duration is a verified value
    (NT, '0_02:00:00\\n1_00:00:00")', '0_02:00:00")', ["namelist_table_current"]),
    # the radiation-replay check (tests/radiation_replay, RG's calls at 00:00 and 00:30): a dump whose cldfrac record
    # holds o3vmr (write_radiation_dump) must fail the replay
    (RDRV, '  reals("cldfrac", levels(cv(m.cldfrac), nV));',
     '  reals("cldfrac", levels(cv(m.o3vmr), nV));', ["radiation_replay_RG_0", "radiation_replay_RG_1800"]),
]

MUTATIONS = [(p, o, n, [t for t in ts if t not in EQUIVALENT_ON_SC.get((p, n), [])]) for p, o, n, ts in MUTATIONS]


def test_status(out, t):
    """The status CTest prints on the named test's own line ("Passed", "***Failed", "***Not Run", "***Exception:
    SegFault", "Subprocess aborted***Exception:", ...), or None if the test has no such line."""
    m = re.search(rf"^\s*\d+/\d+ Test +#\d+: {re.escape(t)} \.+ *(.*?)\s+[\d.]+ sec$", out, re.M)
    return m.group(1).strip() if m else None


def killed(status):
    """A test kills the mutant only if it ran and failed. "Not Run" (a fixture setup the mutant broke, a missing
    executable, a disabled test) is a check that did not run: not a kill (CONSTITUTION.md 3.2). Nor is "Skipped" (the
    test reported that it could not run) or "Timeout" (a stalled machine looks the same as a mutant that hangs)."""
    return (status is not None and status != "Passed" and "Not Run" not in status and "Skipped" not in status
            and "Timeout" not in status)


def is_kill(rc, status, t, broken):
    """CTest failed (rc != 0), the test's own line shows a failure, and the test passed on the unmutated source."""
    return rc != 0 and killed(status) and t not in broken


def not_passing(out, tests):
    """The tests that did not pass on the unmutated source: their failure under a mutant would not show a kill."""
    return {t for t in tests if test_status(out, t) != "Passed"}


def parse_selection(sel, n):
    """The mutant indices of the command-line selection: 'k' (k .. n-1) or 'i,j,...' (two or more distinct indices).
    Anything else (an empty item, a trailing ',', a sign, an index past the list, a repeated index) is refused."""
    if re.fullmatch(r"[0-9]+", sel):
        if int(sel) >= n:
            raise SystemExit(f"mutation_check.py: first mutant {sel} is past the last one ({n - 1})")
        return list(range(int(sel), n))
    if not re.fullmatch(r"[0-9]+(,[0-9]+)+", sel):
        raise SystemExit(f"mutation_check.py: mutant selection {sel!r} is neither an index nor a list i,j,... of "
                         f"indices")
    idx = [int(i) for i in sel.split(",")]
    for i in idx:
        if i >= n:
            raise SystemExit(f"mutation_check.py: mutant {i} is past the last one ({n - 1})")
    if len(set(idx)) != len(idx):
        raise SystemExit(f"mutation_check.py: mutant selection {sel!r} repeats an index")
    return idx


def ctest(build, regex):
    """CTest on the tests the regex selects, one at a time (-j 1 overrides CTEST_PARALLEL_LEVEL)"""
    return subprocess.run(["ctest", "--test-dir", str(build), "-j", "1", "-R", regex, "-FA", "inputs"],
                          capture_output=True, text=True)


def any_of(tests):
    return "^(" + "|".join(re.escape(t) for t in tests) + ")$"


def main():
    if len(sys.argv) not in (3, 4):
        raise SystemExit(__doc__)
    src, build = Path(sys.argv[1]), Path(sys.argv[2])
    chosen = [MUTATIONS[i] for i in parse_selection(sys.argv[3] if len(sys.argv) > 3 else "0", len(MUTATIONS))]
    covered, bad = set(), []
    # every test a mutant must break has to pass first, on the unmutated source
    subprocess.run(["cmake", "--build", str(build), "-j4"], capture_output=True, text=True, check=True)
    needed = sorted({t for _, _, _, ts in chosen for t in ts})
    r = ctest(build, any_of(needed))
    broken = not_passing(r.stdout, needed)
    for t in sorted(broken):
        bad.append(f"{t} does not pass on the unmutated source (status {test_status(r.stdout, t)!r}): no kill counted")
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
            r = ctest(build, any_of(tests)) if tests else None
            for t in tests:
                status = test_status(r.stdout, t)
                failed = is_kill(r.returncode, status, t, broken)
                print(f"{'KILLED' if failed else 'SURVIVED'} {t:22s} by {path}: {new.strip()[:70]}")
                if failed:
                    covered.add(t)
                else:
                    bad.append(f"{t} survived {path}: {new.strip()} (status {status!r})")
        finally:
            p.write_text(orig)
    subprocess.run(["cmake", "--build", str(build), "-j4"], capture_output=True, text=True, check=True)
    # the restored source: a test that killed a mutant must pass again (a stale object, a scratch file a mutant left
    # behind or a change of the machine would otherwise show as kills)
    if covered:
        r = ctest(build, any_of(sorted(covered)))
        for t in sorted(not_passing(r.stdout, covered)):
            bad.append(f"{t} does not pass on the restored source (status {test_status(r.stdout, t)!r}): its kills "
                       f"in this run are not trusted")
    for b in bad:
        print("PROBLEM", b)
    print(f"{len(covered)} tests killed by at least one mutant")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
