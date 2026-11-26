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
    # --- OPM05: every initdump field (atm_mpas_init_block at t = 0) against mpas842-dbl-mw-003 / sc-002 / jw-002
    # atm_compute_solve_diagnostics
    (S, "h_edge(k, iEdge) = 0.5 * (h(k, cell1) + h(k, cell2));", "h_edge(k, iEdge) = 0.5 * (h(k, cell1) + h(k, cell1));",
     D("rho_edge", "MW JW")),
    (S, "ke_edge(k, iEdge) = efac * (u(k, iEdge) * u(k, iEdge));", "ke_edge(k, iEdge) = (efac * u(k, iEdge)) * u(k, iEdge);",
     D("ke", "MW JW")),
    (S, "    for (int i = 1; i <= vertexDegree; ++i) {", "    for (int i = vertexDegree; i >= 1; --i) {",
     D("vorticity pv_vertex", "MW JW")),
    (S, "divergence(k, iCell) = divergence(k, iCell) * r;", "divergence(k, iCell) = divergence(k, iCell) / m.areaCell(iCell);",
     D("divergence", "MW JW")),
    (S, "(ke_edge(k, e1) + ke_edge(k, e2) + ke_edge(k, e3)) * r;", "(ke_edge(k, e3) + ke_edge(k, e2) + ke_edge(k, e1)) * r;",
     D("ke", "MW JW")),
    (S, "(1. - ke_fact) * m.kiteAreasOnVertex(j, iVertex) * ke_vertex(k, iVertex) * r;",
     "(1. - ke_fact) * m.kiteAreasOnVertex(j, iVertex) * (ke_vertex(k, iVertex) * r);", D("ke", "MW JW")),
    (S, "      for (int i = 1; i <= m.nEdgesOnEdge(iEdge); ++i) {", "      for (int i = m.nEdgesOnEdge(iEdge); i >= 1; --i) {",
     D("v", "MW JW")),
    (S, "pv_edge(k, iEdge) = 0.5 * (pv_vertex(k, m.verticesOnEdge(1, iEdge)) + pv_vertex(k, m.verticesOnEdge(2, iEdge)));",
     "pv_edge(k, iEdge) = 0.5 * (pv_vertex(k, m.verticesOnEdge(1, iEdge)) + pv_vertex(k, m.verticesOnEdge(1, iEdge)));",
     D("pv_edge", "MW JW")),
    (S, "pv_cell(k, iCell) + m.kiteAreasOnVertex(j, iVertex) * pv_vertex(k, iVertex) * r;",
     "pv_cell(k, iCell) + m.kiteAreasOnVertex(j, iVertex) * (pv_vertex(k, iVertex) * r);", D("pv_cell gradPVn", "MW JW")),
    (S, "      gradPVt(k, iEdge) = (pv_vertex(k, v2) - pv_vertex(k, v1)) * r1;",
     "      gradPVt(k, iEdge) = (pv_vertex(k, v2) - pv_vertex(k, v1)) * r2;", D("gradPVt", "MW JW")),
    (S, "  const real r = config_apvm_upwinding * dt;", "  const real r = config_apvm_upwinding;", D("pv_edge", "MW JW")),
    (S, "  if (!(config_apvm_upwinding > 0.0))", "  if (false)", ["err_apvm_off"]),
    (S, "  if (!(config_apvm_upwinding > 0.0))", "  if (!(config_apvm_upwinding > 1.0))",
     ["initdump_run_MW", "initdump_run_SC", "initdump_run_JW"]),
    # mpas_reconstruct
    (RC, "uReconstructX(k, iCell) + coeffs_reconstruct(1, i, iCell) * u(k, iEdge);",
     "uReconstructX(k, iCell) + coeffs_reconstruct(2, i, iCell) * u(k, iEdge);", D("uReconstructX uReconstructZonal", "MW JW")),
    (RC, "uReconstructY(k, iCell) + coeffs_reconstruct(2, i, iCell) * u(k, iEdge);",
     "uReconstructY(k, iCell) + coeffs_reconstruct(1, i, iCell) * u(k, iEdge);", D("uReconstructY uReconstructMeridional", "MW JW")),
    (RC, "uReconstructZ(k, iCell) + coeffs_reconstruct(3, i, iCell) * u(k, iEdge);",
     "uReconstructZ(k, iCell) + coeffs_reconstruct(1, i, iCell) * u(k, iEdge);", D("uReconstructZ", "MW JW")),
    (RC, "        uReconstructZonal(k, iCell) = uReconstructX(k, iCell);", "        uReconstructZonal(k, iCell) = uReconstructY(k, iCell);",
     D("uReconstructZonal", "MW SC")),
    (RC, "        uReconstructMeridional(k, iCell) = uReconstructY(k, iCell);",
     "        uReconstructMeridional(k, iCell) = uReconstructX(k, iCell);", D("uReconstructMeridional", "MW SC")),
    (RC, "      libm::sincos(m.latCell(iCell), &slat, &clat);",
     "      slat = libm::sin(m.latCell(iCell)); clat = libm::cos(m.latCell(iCell));", D("uReconstructZonal uReconstructMeridional", "JW")),
    # the arrays OPM03 could not compare
    (A, "m.invAreaCell(iCell) = 1.0 / m.areaCell(iCell);",
     "m.invAreaCell(iCell) = 1.0 / std::sqrt(m.areaCell(iCell)) / std::sqrt(m.areaCell(iCell));", D("invAreaCell", "MW JW")),
    (A, "m.invDvEdge(iEdge) = 1.0 / m.dvEdge(iEdge);",
     "m.invDvEdge(iEdge) = 1.0 / std::sqrt(m.dvEdge(iEdge)) / std::sqrt(m.dvEdge(iEdge));", D("invDvEdge", "MW JW")),
    (A, "m.invDcEdge(iEdge) = 1.0 / m.dcEdge(iEdge);",
     "m.invDcEdge(iEdge) = 1.0 / std::sqrt(m.dcEdge(iEdge)) / std::sqrt(m.dcEdge(iEdge));", D("invDcEdge", "MW JW")),
    (A, "m.invAreaTriangle(iVertex) = 1.0 / m.areaTriangle(iVertex);",
     "m.invAreaTriangle(iVertex) = 1.0 / std::sqrt(m.areaTriangle(iVertex)) / std::sqrt(m.areaTriangle(iVertex));",
     D("invAreaTriangle", "MW JW")),
    (A, "        m.edgesOnCell_sign(i, iCell) = side == 1 ? 1.0 : -1.0;", "        m.edgesOnCell_sign(i, iCell) = 1.0;",
     D("edgesOnCell_sign", "MW JW")),
    (A, "          m.edgesOnVertex_sign(i, iVtx) = -1.0;", "          m.edgesOnVertex_sign(i, iVtx) = 1.0;",
     D("edgesOnVertex_sign", "MW JW")),
    (A, "          m.zb_cell(k, i, iCell) = m.zb(k, side, e);", "          m.zb_cell(k, i, iCell) = m.zb(k, 1, e);",
     D("zb_cell", "MW JW")),
    (A, "          m.zb3_cell(k, i, iCell) = m.zb3(k, side, e);", "          m.zb3_cell(k, i, iCell) = m.zb3(k, 1, e);",
     D("zb3_cell", "MW JW")),
    (A, "            m.kiteForCell(i, iCell) = j;", "            m.kiteForCell(i, iCell) = 1;", D("kiteForCell", "MW JW")),
    (A, "    for (int iCell = 1; iCell <= n; ++iCell) m.advCellsForEdge(iCell, iEdge) = cell_list[iCell];",
     "    for (int iCell = 1; iCell <= n; ++iCell) m.advCellsForEdge(iCell, iEdge) = cell_list[1];", D("advCellsForEdge", "MW JW")),
    (A, "    m.nAdvCellsForEdge(iEdge) = n;", "    m.nAdvCellsForEdge(iEdge) = n - 1;", D("nAdvCellsForEdge", "MW JW")),
    (A, "    j_in = slot(cell1);  // :1513-1517\n    adv(j_in, iEdge) = adv(j_in, iEdge) + 0.5;",
     "    j_in = slot(cell1);  // :1513-1517\n    adv(j_in, iEdge) = adv(j_in, iEdge) + 0.25;", D("adv_coefs", "MW JW")),
    (A, "m.adv_coefs_3rd.data()[k] = c * m.adv_coefs_3rd.data()[k];", "m.adv_coefs_3rd.data()[k] = m.adv_coefs_3rd.data()[k];",
     D("adv_coefs_3rd", "MW JW")),
    (A, "      if (m.cellsOnEdge(1, m.edgesOnCell(i, iCell)) != iCell) {", "      if (false) {",
     D("deformation_coef_c2 deformation_coef_s2 deformation_coef_cs", "MW JW")),
    (A, "      thetat[i] = libm::atan2(dy, dx) - pii / 2.;", "      thetat[i] = libm::atan2(dy, dx) + pii / 2.;",
     D("deformation_coef_c deformation_coef_s", "MW JW")),
    (A, "  fill_all(m.meshScalingDel2, 1.0);", "  fill_all(m.meshScalingDel2, 0.5);", D("meshScalingDel2", "MW JW")),
    (A, "  fill_all(m.meshScalingDel4, 1.0);", "  fill_all(m.meshScalingDel4, 0.5);", D("meshScalingDel4", "MW JW")),
    (A, "  fill_all(m.meshScalingRegionalCell, 1.0);", "  fill_all(m.meshScalingRegionalCell, 0.5);",
     D("meshScalingRegionalCell", "MW JW")),
    (A, "  fill_all(m.meshScalingRegionalEdge, 1.0);", "  fill_all(m.meshScalingRegionalEdge, 0.5);",
     D("meshScalingRegionalEdge", "MW JW")),
    (A, "      epssm_coeff_u[k] = min_coeff;", "      epssm_coeff_u[k] = max_coeff;", D("etp etm", "MW JW")),
    (A, "      epssm_coeff_w[k] = min_coeff;", "      epssm_coeff_w[k] = max_coeff;", D("ewp ewm", "MW JW")),
    (A, "    if (m.bdyMaskCell.data()[k] > nRelaxZone)", "    if (m.bdyMaskCell.data()[k] >= 0)", D("specZoneMaskCell", "MW JW")),
    (A, "    if (m.bdyMaskEdge.data()[k] > nRelaxZone)", "    if (m.bdyMaskEdge.data()[k] >= 0)", D("specZoneMaskEdge", "MW JW")),
    (A, "    if (m.bdyMaskVertex.data()[k] > nRelaxZone)", "    if (m.bdyMaskVertex.data()[k] >= 0)", D("specZoneMaskVertex", "MW JW")),
    (A, "  fill_all(m.nearestRelaxationCell, nCells + 1);", "  fill_all(m.nearestRelaxationCell, nCells);",
     D("nearestRelaxationCell", "MW JW")),
    (A, "      theta_m(k, iCell) = s.theta(k, iCell) * (1. + rvord * s.scalars(index_qv, k, iCell));",
     "      theta_m(k, iCell) = s.theta_base(k, iCell) * (1. + rvord * s.scalars(index_qv, k, iCell));", D("theta_m rtheta_p", "MW JW")),
    (A, "      rho_zz(k, iCell) = s.rho(k, iCell) / zz(k, iCell);", "      rho_zz(k, iCell) = s.rho(k, iCell) * (1.0 / zz(k, iCell));",
     D("rho_zz ru rho_p", "MW JW")),
    (A, "      d.ru(k, iEdge) = 0.5 * s.u(k, iEdge) * (rho_zz(k, cell1) + rho_zz(k, cell2));",
     "      d.ru(k, iEdge) = 0.5 * (s.u(k, iEdge) * (rho_zz(k, cell1) + rho_zz(k, cell2)));", D("ru", "MW JW")),
    (A, "(m.zb_cell(k, i, iCell) + std::copysign(1.0, flux) * m.zb3_cell(k, i, iCell))",
     "(m.zb_cell(k, i, iCell) + m.zb3_cell(k, i, iCell))", D("rw", "MW JW")),
    (A, "d.rho_p(k, iCell) = rho_zz(k, iCell) - s.rho_base(k, iCell);", "d.rho_p(k, iCell) = s.rho(k, iCell) - s.rho_base(k, iCell);",
     D("rho_p", "MW JW")),
    (A, "      d.exner(k, iCell) = vector ? libm::mvec_pow(xe, rcv) : libm::pow(xe, rcv);",
     "      d.exner(k, iCell) = libm::pow(xe, rcv);", D("exner pressure_p", "MW JW")),
    (A, "      d.exner_base(k, iCell) = vector ? libm::mvec_pow(xb, rcv) : libm::pow(xb, rcv);",
     "      d.exner_base(k, iCell) = libm::pow(xb, rcv);", D("exner_base pressure_base", "MW JW")),
    (A, "    for (int k = 1; k <= nVertLevels; ++k) d.rtheta_base(k, iCell) = s.theta_base(k, iCell) * s.rho_base(k, iCell);",
     "    for (int k = 1; k <= nVertLevels; ++k) d.rtheta_base(k, iCell) = s.rho_base(k, iCell) * s.theta_base(k, iCell) * 1.0000000000000002;",
     D("rtheta_base", "MW JW")),
    (W, "    std::vector<char> v(ref.c->begin(), ref.c->end());", "    std::vector<char> v(ref.c->rbegin(), ref.c->rend());",
     D("xtime", "MW JW")),
    ("tests/CMakeLists.txt", "  xtime\n  invAreaCell", "  invAreaCell", ["initdump_vars_MW", "initdump_vars_SC", "initdump_vars_JW"]),
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
