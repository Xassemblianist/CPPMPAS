// SPDX-License-Identifier: GPL-3.0-or-later
#include "core_atmosphere/atm_init.hpp"

#include <algorithm>
#include <cmath>

#include "core/error.hpp"
#include "core/libm.hpp"
#include "operators/geometry_utils.hpp"
#include "operators/vector_operations.hpp"

namespace cppmpas {

namespace {

template <class T> void fill_all(FArray<T>& a, T v) { std::fill(a.data(), a.data() + a.size(), v); }

// Physical constants, src/framework/mpas_constants.F:46-52 (gfortran folds the parameter expressions with one
// correctly rounded operation each, as the C++ compiler does here)
constexpr real rgas = 287.0;
constexpr real rv = 461.6;
constexpr real cp = 7.0 * rgas / 2.0;
constexpr real rvord = rv / rgas;

// pi as acos(m1) (CORE:1292-1293) and 2.*asin(1.0) (CORE:1698): folded at compile time in the reference (no acos or
// asin call in atm_compute_damping_coefs; atm_initialize_deformation_weights calls asin only in the geometry
// functions) to the double nearest pi
constexpr real kPi = 0x1.921fb54442d18p+1;

}  // namespace

void allocate_init_fields(Mesh& m, State& s, Diag& d) {
  if (m.nVertLevels <= 0) throw Error("atm_mpas_init_block: the mesh has no vertical levels (not an init file)");
  const int nC = m.nCells + 1, nE = m.nEdges + 1, nV = m.nVertices + 1, nL = m.nVertLevels, nL1 = m.nVertLevels + 1;
  // Registry default value 0 / 0.0 for all of them (no default_value, gen_inc.c:615-630; specZoneMask* "0.0")
  m.invAreaCell.allocate({nC}, 0.0);
  m.invDvEdge.allocate({nE}, 0.0);
  m.invDcEdge.allocate({nE}, 0.0);
  m.invAreaTriangle.allocate({nV}, 0.0);
  m.edgesOnCell_sign.allocate({m.maxEdges, nC}, 0.0);
  m.edgesOnVertex_sign.allocate({m.vertexDegree, nV}, 0.0);
  m.zb_cell.allocate({nL1, m.maxEdges, nC}, 0.0);
  m.zb3_cell.allocate({nL1, m.maxEdges, nC}, 0.0);
  m.kiteForCell.allocate({m.maxEdges, nC}, 0);
  m.advCellsForEdge.allocate({15, nE}, 0);
  m.nAdvCellsForEdge.allocate({nE}, 0);
  m.adv_coefs.allocate({15, nE}, 0.0);
  m.adv_coefs_3rd.allocate({15, nE}, 0.0);
  for (auto* a : {&m.deformation_coef_c2, &m.deformation_coef_s2, &m.deformation_coef_cs, &m.deformation_coef_c,
                  &m.deformation_coef_s})
    a->allocate({m.maxEdges, nC}, 0.0);
  m.meshScalingDel2.allocate({nE}, 0.0);
  m.meshScalingDel4.allocate({nE}, 0.0);
  m.meshScalingRegionalCell.allocate({nC}, 0.0);
  m.meshScalingRegionalEdge.allocate({nE}, 0.0);
  m.etp.allocate({nL}, 0.0);
  m.etm.allocate({nL}, 0.0);
  m.ewp.allocate({nL1}, 0.0);
  m.ewm.allocate({nL1}, 0.0);
  m.specZoneMaskCell.allocate({nC}, 0.0);
  m.specZoneMaskEdge.allocate({nE}, 0.0);
  m.specZoneMaskVertex.allocate({nV}, 0.0);
  m.nearestRelaxationCell.allocate({nC}, 0);
  s.theta_m.allocate({nL, nC}, 0.0);
  s.rho_zz.allocate({nL, nC}, 0.0);
  d.ru.allocate({nL, nE}, 0.0);
  d.rw.allocate({nL1, nC}, 0.0);
  for (auto* a : {&d.rho_p, &d.rtheta_base, &d.rtheta_p, &d.exner, &d.exner_base, &d.pressure_p, &d.pressure_base})
    a->allocate({nL, nC}, 0.0);
}

// mpas_atm_core.F:1184-1271
void atm_compute_signs(Mesh& m, std::vector<std::string>& log) {
  const int nCells = m.nCells, nVertices = m.nVertices, nEdges = m.nEdges, vertexDegree = m.vertexDegree;
  const int nL1 = m.nVertLevels + 1;
  for (int iVtx = 1; iVtx <= nVertices; ++iVtx) {  // :1218-1230
    for (int i = 1; i <= vertexDegree; ++i) {
      if (m.edgesOnVertex(i, iVtx) <= nEdges) {
        if (iVtx == m.verticesOnEdge(2, m.edgesOnVertex(i, iVtx)))
          m.edgesOnVertex_sign(i, iVtx) = 1.0;
        else
          m.edgesOnVertex_sign(i, iVtx) = -1.0;
      } else {
        m.edgesOnVertex_sign(i, iVtx) = 0.0;
      }
    }
  }
  for (int iCell = 1; iCell <= nCells; ++iCell) {  // :1232-1248
    for (int i = 1; i <= m.nEdgesOnCell(iCell); ++i) {
      const int e = m.edgesOnCell(i, iCell);
      if (e <= nEdges) {
        const int side = (iCell == m.cellsOnEdge(1, e)) ? 1 : 2;
        m.edgesOnCell_sign(i, iCell) = side == 1 ? 1.0 : -1.0;
        for (int k = 1; k <= nL1; ++k) {
          m.zb_cell(k, i, iCell) = m.zb(k, side, e);
          m.zb3_cell(k, i, iCell) = m.zb3(k, side, e);
        }
      } else {
        m.edgesOnCell_sign(i, iCell) = 0.0;
      }
    }
  }
  for (int iCell = 1; iCell <= nCells; ++iCell) {  // :1250-1268
    for (int i = 1; i <= m.nEdgesOnCell(iCell); ++i) {
      const int iVtx = m.verticesOnCell(i, iCell);
      if (iVtx <= nVertices) {
        int j = 1;
        for (; j <= vertexDegree; ++j) {
          if (iCell == m.cellsOnVertex(j, iVtx)) {
            m.kiteForCell(i, iCell) = j;
            break;
          }
        }
        if (j > vertexDegree) log.push_back("ERROR: Unexpected error while identifying kiteForCell");  // :1262-1264
      } else {
        m.kiteForCell(i, iCell) = 1;
      }
    }
  }
}

// mpas_atm_core.F:459-473
void atm_compute_inverses(Mesh& m) {
  for (int iCell = 1; iCell <= m.nCells; ++iCell) m.invAreaCell(iCell) = 1.0 / m.areaCell(iCell);
  for (int iEdge = 1; iEdge <= m.nEdges; ++iEdge) m.invDvEdge(iEdge) = 1.0 / m.dvEdge(iEdge);
  for (int iEdge = 1; iEdge <= m.nEdges; ++iEdge) m.invDcEdge(iEdge) = 1.0 / m.dcEdge(iEdge);
  for (int iVertex = 1; iVertex <= m.nVertices; ++iVertex)
    m.invAreaTriangle(iVertex) = 1.0 / m.areaTriangle(iVertex);
}

// mpas_atm_core.F:1620-1850
void atm_initialize_deformation_weights(Mesh& m, int nCells, bool on_a_sphere, real sphere_radius) {
  const real x_period = m.x_period, y_period = m.y_period;  // :1664-1665
  fill_all(m.deformation_coef_c2, 0.0);                     // :1692-1696
  fill_all(m.deformation_coef_s2, 0.0);
  fill_all(m.deformation_coef_cs, 0.0);
  fill_all(m.deformation_coef_c, 0.0);
  fill_all(m.deformation_coef_s, 0.0);
  const real pii = kPi;  // :1698
  std::vector<real> theta_abs(static_cast<std::size_t>(nCells) + 1);
  real xc[26], yc[26], zc[26], thetav[26], thetat[26], dl_sphere[26], xp[26], yp[26];  // dimension(25), 1-based
  int cell_list[26];

  for (int iCell = 1; iCell <= nCells; ++iCell) {
    cell_list[1] = iCell;  // :1702-1706
    for (int i = 2; i <= m.nEdgesOnCell(iCell) + 1; ++i) cell_list[i] = m.cellsOnCell(i - 1, iCell);
    const int n = m.nEdgesOnCell(iCell) + 1;
    if (n > 25) throw Error("atm_initialize_deformation_weights: nEdgesOnCell > 24 overflows MPAS's dimension(25)");
    bool do_the_cell = true;  // :1710-1716
    for (int i = 1; i <= n; ++i)
      if (cell_list[i] > nCells) do_the_cell = false;
    if (!do_the_cell) continue;

    if (on_a_sphere) {
      xc[1] = m.xCell(iCell) / sphere_radius;  // :1725-1734
      yc[1] = m.yCell(iCell) / sphere_radius;
      zc[1] = m.zCell(iCell) / sphere_radius;
      for (int i = 2; i <= n; ++i) {
        const int iv = m.verticesOnCell(i - 1, iCell);
        xc[i] = m.xVertex(iv) / sphere_radius;
        yc[i] = m.yVertex(iv) / sphere_radius;
        zc[i] = m.zVertex(iv) / sphere_radius;
      }
      if (zc[1] == 1.0) {  // :1742-1750
        theta_abs[iCell] = pii / 2.;
      } else {
        theta_abs[iCell] = pii / 2. - mpas_sphere_angle(xc[1], yc[1], zc[1], xc[2], yc[2], zc[2], 0.0, 0.0, 1.0);
      }
      thetat[1] = theta_abs[iCell];  // :1760-1772
      for (int i = 1; i <= n - 1; ++i) {
        int ip2 = i + 2;
        if (ip2 > n) ip2 = 2;
        thetav[i] = mpas_sphere_angle(xc[1], yc[1], zc[1], xc[i + 1], yc[i + 1], zc[i + 1], xc[ip2], yc[ip2], zc[ip2]);
        dl_sphere[i] = sphere_radius * mpas_arc_length(xc[1], yc[1], zc[1], xc[i + 1], yc[i + 1], zc[i + 1]);
        if (i > 1) thetat[i] = thetat[i - 1] + thetav[i - 1];
      }
      // :1776-1779. gfortran vectorised this loop: pairs i = (1,2), (3,4), ... call _ZGVbN2v_cos and _ZGVbN2v_sin
      // (one lane each); an odd last i calls sincos@plt (reference executable, atm_initialize_deformation_weights
      // +0x1a90 ... +0x1b23).
      const int cnt = n - 1;
      int i = 1;
      for (; i + 1 <= cnt; i += 2) {
        for (int l = i; l <= i + 1; ++l) {
          xp[l] = libm::mvec_cos(thetat[l]) * dl_sphere[l];
          yp[l] = libm::mvec_sin(thetat[l]) * dl_sphere[l];
        }
      }
      if (i <= cnt) {
        real s, c;
        libm::sincos(thetat[i], &s, &c);
        xp[i] = c * dl_sphere[i];
        yp[i] = s * dl_sphere[i];
      }
    } else {
      for (int i = 1; i <= n - 1; ++i) {  // :1785-1789
        const int iv = m.verticesOnCell(i, iCell);
        xp[i] = mpas_fix_periodicity(m.xVertex(iv), m.xCell(iCell), x_period) - m.xCell(iCell);
        yp[i] = mpas_fix_periodicity(m.yVertex(iv), m.yCell(iCell), y_period) - m.yCell(iCell);
      }
      for (int i = 1; i <= n - 1; ++i) {  // :1791-1796, atan2@plt
        const int ie = m.edgesOnCell(i, iCell);
        const real xe = mpas_fix_periodicity(m.xEdge(ie), m.xCell(iCell), x_period) - m.xCell(iCell);
        const real ye = mpas_fix_periodicity(m.yEdge(ie), m.yCell(iCell), y_period) - m.yCell(iCell);
        thetat[i] = libm::atan2(ye, xe);
      }
      theta_abs[iCell] = thetat[1];  // :1798
    }

    thetat[1] = theta_abs[iCell];  // :1804-1813 (values overwritten at :1821 before any read)
    for (int i = 2; i <= n - 1; ++i) {
      int ip1 = i + 1;
      if (ip1 == n) ip1 = 1;
      thetat[i] = mpas_plane_angle(0.0, 0.0, 0.0, xp[i] - xp[i - 1], yp[i] - yp[i - 1], 0.0, xp[ip1] - xp[i],
                                   yp[ip1] - yp[i], 0.0, 0.0, 0.0, 1.0);
      thetat[i] = thetat[i] + thetat[i - 1];
    }

    real area_cell = 0.;  // :1815-1823, atan2@plt
    for (int i = 1; i <= n - 1; ++i) {
      int ip1 = i + 1;
      if (ip1 == n) ip1 = 1;
      const real dx = xp[ip1] - xp[i];
      const real dy = yp[ip1] - yp[i];
      area_cell = area_cell + 0.25 * (xp[i] + xp[ip1]) * (yp[ip1] - yp[i]) - 0.25 * (yp[i] + yp[ip1]) * (xp[ip1] - xp[i]);
      thetat[i] = libm::atan2(dy, dx) - pii / 2.;
    }

    for (int i = 1; i <= n - 1; ++i) {  // :1827-1847; sin / cos of thetat(i) fused into one sincos@plt
      int ip1 = i + 1;
      if (ip1 == n) ip1 = 1;
      const real dl = std::sqrt((xp[ip1] - xp[i]) * (xp[ip1] - xp[i]) + (yp[ip1] - yp[i]) * (yp[ip1] - yp[i]));
      real s, c;
      libm::sincos(thetat[i], &s, &c);
      const real sint2 = s * s;
      const real cost2 = c * c;
      const real sint_cost = s * c;
      m.deformation_coef_c2(i, iCell) = dl * cost2 / area_cell;
      m.deformation_coef_s2(i, iCell) = dl * sint2 / area_cell;
      m.deformation_coef_cs(i, iCell) = dl * sint_cost / area_cell;
      m.deformation_coef_c(i, iCell) = dl * c / area_cell;
      m.deformation_coef_s(i, iCell) = dl * s / area_cell;
      if (m.cellsOnEdge(1, m.edgesOnCell(i, iCell)) != iCell) {
        m.deformation_coef_c2(i, iCell) = -m.deformation_coef_c2(i, iCell);
        m.deformation_coef_s2(i, iCell) = -m.deformation_coef_s2(i, iCell);
        m.deformation_coef_cs(i, iCell) = -m.deformation_coef_cs(i, iCell);
      }
    }
  }
}

// mpas_atm_core.F:1391-1536
void atm_adv_coef_compression(Mesh& m) {
  const int nCells = m.nCells, nEdges = m.nEdges;
  auto& adv = m.adv_coefs;
  auto& adv3 = m.adv_coefs_3rd;
  int cell_list[21];
  for (int iEdge = 1; iEdge <= nEdges; ++iEdge) {
    m.nAdvCellsForEdge(iEdge) = 0;
    const int cell1 = m.cellsOnEdge(1, iEdge);
    const int cell2 = m.cellsOnEdge(2, iEdge);
    if (!(cell1 <= nCells || cell2 <= nCells)) continue;  // :1431
    cell_list[1] = cell1;
    cell_list[2] = cell2;
    int n = 2;
    for (int i = 1; i <= m.nEdgesOnCell(cell1); ++i) {  // :1438-1443
      if (m.cellsOnCell(i, cell1) != cell2) {
        n = n + 1;
        if (n > 20) throw Error("atm_adv_coef_compression: more than 20 cells overflow MPAS's cell_list(20)");
        cell_list[n] = m.cellsOnCell(i, cell1);
      }
    }
    for (int iCell = 1; iCell <= m.nEdgesOnCell(cell2); ++iCell) {  // :1446-1456
      bool addcell = true;
      for (int i = 1; i <= n; ++i)
        if (cell_list[i] == m.cellsOnCell(iCell, cell2)) addcell = false;
      if (addcell) {
        n = n + 1;
        if (n > 20) throw Error("atm_adv_coef_compression: more than 20 cells overflow MPAS's cell_list(20)");
        cell_list[n] = m.cellsOnCell(iCell, cell2);
      }
    }
    if (n > 15) throw Error("atm_adv_coef_compression: more than 15 cells overflow advCellsForEdge(FIFTEEN, nEdges)");
    m.nAdvCellsForEdge(iEdge) = n;  // :1459-1462
    for (int iCell = 1; iCell <= n; ++iCell) m.advCellsForEdge(iCell, iEdge) = cell_list[iCell];

    for (int j = 1; j <= 15; ++j) {  // :1466-1467
      adv(j, iEdge) = 0.;
      adv3(j, iEdge) = 0.;
    }
    auto slot = [&](int c) {  // the last j with cell_list(j) == c (0: none)
      int j_in = 0;
      for (int j = 1; j <= n; ++j)
        if (cell_list[j] == c) j_in = j;
      if (j_in == 0) throw Error("atm_adv_coef_compression: cell not in cell_list (MPAS would index slot 0)");
      return j_in;
    };
    int j_in = slot(cell1);  // :1472-1477
    adv(j_in, iEdge) = adv(j_in, iEdge) + m.deriv_two(1, 1, iEdge);
    adv3(j_in, iEdge) = adv3(j_in, iEdge) + m.deriv_two(1, 1, iEdge);
    for (int iCell = 1; iCell <= m.nEdgesOnCell(cell1); ++iCell) {  // :1479-1486
      j_in = slot(m.cellsOnCell(iCell, cell1));
      adv(j_in, iEdge) = adv(j_in, iEdge) + m.deriv_two(iCell + 1, 1, iEdge);
      adv3(j_in, iEdge) = adv3(j_in, iEdge) + m.deriv_two(iCell + 1, 1, iEdge);
    }
    j_in = slot(cell2);  // :1491-1496
    adv(j_in, iEdge) = adv(j_in, iEdge) + m.deriv_two(1, 2, iEdge);
    adv3(j_in, iEdge) = adv3(j_in, iEdge) - m.deriv_two(1, 2, iEdge);
    for (int iCell = 1; iCell <= m.nEdgesOnCell(cell2); ++iCell) {  // :1498-1505
      j_in = slot(m.cellsOnCell(iCell, cell2));
      adv(j_in, iEdge) = adv(j_in, iEdge) + m.deriv_two(iCell + 1, 2, iEdge);
      adv3(j_in, iEdge) = adv3(j_in, iEdge) - m.deriv_two(iCell + 1, 2, iEdge);
    }
    for (int j = 1; j <= n; ++j) {  // :1506-1509
      adv(j, iEdge) = -((m.dcEdge(iEdge) * m.dcEdge(iEdge)) * adv(j, iEdge) / 12.);
      adv3(j, iEdge) = -((m.dcEdge(iEdge) * m.dcEdge(iEdge)) * adv3(j, iEdge) / 12.);
    }
    j_in = slot(cell1);  // :1513-1517
    adv(j_in, iEdge) = adv(j_in, iEdge) + 0.5;
    j_in = slot(cell2);  // :1519-1523
    adv(j_in, iEdge) = adv(j_in, iEdge) + 0.5;
    for (int j = 1; j <= n; ++j) {  // :1527-1530
      adv(j, iEdge) = m.dvEdge(iEdge) * adv(j, iEdge);
      adv3(j, iEdge) = m.dvEdge(iEdge) * adv3(j, iEdge);
    }
  }
}

// mpas_atm_core.F:1539-1558: whole arrays, garbage element and unused slots included
void atm_couple_coef_3rd_order(Mesh& m, const InitConfig& cfg) {
  const real c = cfg.config_coef_3rd_order;
  for (std::size_t k = 0; k < m.adv_coefs_3rd.size(); ++k) m.adv_coefs_3rd.data()[k] = c * m.adv_coefs_3rd.data()[k];
  for (std::size_t k = 0; k < m.zb3_cell.size(); ++k) m.zb3_cell.data()[k] = c * m.zb3_cell.data()[k];
}

// mpas_atm_core.F:1124-1181. (a+b)/2.0 is (a+b)*0.5 in the reference (exact). The edge loops call pow@plt; the
// meshScalingRegionalCell loop is vectorised: cells (1,2), (3,4), ... call _ZGVbN2vv_pow (one lane each), an odd
// last cell pow@plt (reference executable, atm_compute_mesh_scaling +0x8e0 ... +0x9ad).
void atm_compute_mesh_scaling(Mesh& m, const InitConfig& cfg) {
  const int nEdges = m.nEdges, nCells = m.nCells;
  fill_all(m.meshScalingDel2, 1.0);  // :1152-1153
  fill_all(m.meshScalingDel4, 1.0);
  if (cfg.config_h_ScaleWithMesh) {
    for (int iEdge = 1; iEdge <= nEdges; ++iEdge) {  // :1155-1160
      const int cell1 = m.cellsOnEdge(1, iEdge);
      const int cell2 = m.cellsOnEdge(2, iEdge);
      m.meshScalingDel2(iEdge) = 1.0 / libm::pow((m.meshDensity(cell1) + m.meshDensity(cell2)) / 2.0, 0.25);
      m.meshScalingDel4(iEdge) = 1.0 / libm::pow((m.meshDensity(cell1) + m.meshDensity(cell2)) / 2.0, 0.75);
    }
  }
  fill_all(m.meshScalingRegionalCell, 1.0);  // :1166-1167
  fill_all(m.meshScalingRegionalEdge, 1.0);
  if (cfg.config_h_ScaleWithMesh) {
    for (int iEdge = 1; iEdge <= nEdges; ++iEdge) {  // :1169-1173
      const int cell1 = m.cellsOnEdge(1, iEdge);
      const int cell2 = m.cellsOnEdge(2, iEdge);
      m.meshScalingRegionalEdge(iEdge) = 1.0 / libm::pow((m.meshDensity(cell1) + m.meshDensity(cell2)) / 2.0, 0.25);
    }
    int iCell = 1;  // :1175-1177
    for (; iCell + 1 <= nCells; iCell += 2) {
      m.meshScalingRegionalCell(iCell) = 1.0 / libm::mvec_pow(m.meshDensity(iCell), 0.25);
      m.meshScalingRegionalCell(iCell + 1) = 1.0 / libm::mvec_pow(m.meshDensity(iCell + 1), 0.25);
    }
    if (iCell <= nCells) m.meshScalingRegionalCell(iCell) = 1.0 / libm::pow(m.meshDensity(iCell), 0.25);
  }
}

// mpas_atm_core.F:1274-1388. sin(...)**2.0 and sin(...)**2 are s*s (no pow call in the reference); sin@plt.
void atm_compute_damping_coefs(Mesh& m, const InitConfig& cfg) {
  const int nCells = m.nCells, nVertLevels = m.nVertLevels;
  const real pii = kPi;  // :1292-1293
  const real config_zd = cfg.config_zd, config_xnutr = cfg.config_xnutr;
  fill_all(m.dss, 0.0);  // :1304
  for (int iCell = 1; iCell <= nCells; ++iCell) {  // :1305-1313
    const real zt = m.zgrid(nVertLevels + 1, iCell);
    for (int k = 1; k <= nVertLevels; ++k) {
      const real z = 0.5 * (m.zgrid(k, iCell) + m.zgrid(k + 1, iCell));
      if (z > config_zd) {
        const real s = libm::sin(0.5 * pii * (z - config_zd) / (zt - config_zd));
        m.dss(k, iCell) = config_xnutr * (s * s);
      }
    }
  }
  const real min_coeff = cfg.config_epssm_minimum, max_coeff = cfg.config_epssm_maximum;
  const real transition_lower_bound = cfg.config_epssm_transition_bottom_z;
  const real transition_upper_bound = cfg.config_epssm_transition_top_z;
  std::vector<real> height_u_levels(nVertLevels + 1), epssm_coeff_u(nVertLevels + 1),
      height_w_levels(nVertLevels + 2), epssm_coeff_w(nVertLevels + 2);
  const real transition_width = transition_upper_bound - transition_lower_bound;  // :1332
  std::fill(height_w_levels.begin(), height_w_levels.end(), 0.0);                 // :1337
  for (int k = 1; k <= nVertLevels; ++k) {                                        // :1338-1341
    height_w_levels[k + 1] = height_w_levels[k] + 1.0 / m.rdzw(k);
    height_u_levels[k] = 0.5 * (height_w_levels[k + 1] + height_w_levels[k]);
  }
  for (int k = 1; k <= nVertLevels; ++k) {  // :1354-1368
    if (height_u_levels[k] <= transition_lower_bound) {
      epssm_coeff_u[k] = min_coeff;
    } else if (height_u_levels[k] >= transition_upper_bound) {
      epssm_coeff_u[k] = max_coeff;
    } else {
      const real z = (height_u_levels[k] - transition_lower_bound) / transition_width;
      const real s = libm::sin(0.5 * pii * z);
      epssm_coeff_u[k] = min_coeff + s * s * (max_coeff - min_coeff);
    }
    m.etp(k) = 0.5 * (1.0 + epssm_coeff_u[k]);
    m.etm(k) = 0.5 * (1.0 - epssm_coeff_u[k]);
  }
  for (int k = 1; k <= nVertLevels + 1; ++k) {  // :1369-1383
    if (height_w_levels[k] <= transition_lower_bound) {
      epssm_coeff_w[k] = min_coeff;
    } else if (height_w_levels[k] >= transition_upper_bound) {
      epssm_coeff_w[k] = max_coeff;
    } else {
      const real z = (height_w_levels[k] - transition_lower_bound) / transition_width;
      const real s = libm::sin(0.5 * pii * z);
      epssm_coeff_w[k] = min_coeff + s * s * (max_coeff - min_coeff);
    }
    m.ewp(k) = 0.5 * (1.0 + epssm_coeff_w[k]);
    m.ewm(k) = 0.5 * (1.0 - epssm_coeff_w[k]);
  }
}

// src/core_atmosphere/dynamics/mpas_atm_boundaries.F:659-753 (nRelaxZone = 5, :37)
void mpas_atm_setup_bdy_masks(Mesh& m) {
  constexpr int nRelaxZone = 5;
  const int nCells = m.nCells;
  for (std::size_t k = 0; k < m.bdyMaskCell.size(); ++k)  // :697-699, whole arrays
    if (m.bdyMaskCell.data()[k] > nRelaxZone) m.specZoneMaskCell.data()[k] = 1.0;
  for (std::size_t k = 0; k < m.bdyMaskEdge.size(); ++k)
    if (m.bdyMaskEdge.data()[k] > nRelaxZone) m.specZoneMaskEdge.data()[k] = 1.0;
  for (std::size_t k = 0; k < m.bdyMaskVertex.size(); ++k)
    if (m.bdyMaskVertex.data()[k] > nRelaxZone) m.specZoneMaskVertex.data()[k] = 1.0;
  fill_all(m.nearestRelaxationCell, nCells + 1);  // :702
  auto dist2 = [&](int a, int b) {
    const real dx = m.xCell(a) - m.xCell(b), dy = m.yCell(a) - m.yCell(b), dz = m.zCell(a) - m.zCell(b);
    return dx * dx + dy * dy + dz * dz;
  };
  for (int iCell = 1; iCell <= nCells; ++iCell) {  // :708-722
    if (m.bdyMaskCell(iCell) == nRelaxZone + 1) {
      real dmin = 1.0e36;
      for (int j = 1; j <= m.nEdgesOnCell(iCell); ++j) {
        const int i = m.cellsOnCell(j, iCell);
        if (m.bdyMaskCell(i) == nRelaxZone) {
          const real d = dist2(i, iCell);
          if (d < dmin) {
            dmin = d;
            m.nearestRelaxationCell(iCell) = i;
          }
        }
      }
    }
  }
  for (int iCell = 1; iCell <= nCells; ++iCell) {  // :728-751
    if (m.bdyMaskCell(iCell) == nRelaxZone + 2) {
      real dmin = 1.0e36;
      for (int j = 1; j <= m.nEdgesOnCell(iCell); ++j) {
        const int i = m.cellsOnCell(j, iCell);
        if (m.bdyMaskCell(i) == nRelaxZone + 1) {
          for (int jj = 1; jj <= m.nEdgesOnCell(i); ++jj) {
            const int ii = m.cellsOnCell(jj, i);
            if (m.bdyMaskCell(ii) == nRelaxZone) {
              const real d = dist2(ii, iCell);
              if (d < dmin) {
                dmin = d;
                m.nearestRelaxationCell(iCell) = ii;
              }
            }
          }
        }
      }
    }
  }
}

// src/core_atmosphere/dynamics/mpas_atm_time_integration.F:7440-7665, cellStart..cellEnd = 1..nCells,
// edgeStart..edgeEnd = 1..nEdges (one thread). The exner / exner_base loop (:7633-7638) is vectorised over k:
// levels (1,2), (3,4), ... call _ZGVbN2vv_pow (one lane each), an odd last level pow@plt (reference executable,
// atm_init_coupled_diagnostics +0x41b0 ... +0x43af; the scalar loop +0x4a37 is taken only for nVertLevels = 1 or
// overlapping arrays).
void atm_init_coupled_diagnostics(const Mesh& m, State& s, Diag& d) {
  const int nVertLevels = m.nVertLevels;
  int index_qv = 0;
  for (std::size_t k = 0; k < s.scalar_names.size(); ++k)
    if (s.scalar_names[k] == "qv") index_qv = static_cast<int>(k) + 1;
  if (index_qv == 0) throw Error("atm_init_coupled_diagnostics: qv is not an active scalar");
  constexpr real rcv = rgas / (cp - rgas);  // :7538
  constexpr real p0 = 1.e5;                 // :7539
  constexpr real rgas_p0 = rgas / p0;       // (rgas/p0) of :7635-7636, folded
  const auto& zz = m.zz;
  auto& theta_m = s.theta_m;
  auto& rho_zz = s.rho_zz;

  for (int iCell = 1; iCell <= m.nCells; ++iCell) {  // :7543-7549
    for (int k = 1; k <= nVertLevels; ++k) {
      theta_m(k, iCell) = s.theta(k, iCell) * (1. + rvord * s.scalars(index_qv, k, iCell));
      rho_zz(k, iCell) = s.rho(k, iCell) / zz(k, iCell);
    }
  }
  for (int iEdge = 1; iEdge <= m.nEdges; ++iEdge) {  // :7556-7563
    const int cell1 = m.cellsOnEdge(1, iEdge);
    const int cell2 = m.cellsOnEdge(2, iEdge);
    for (int k = 1; k <= nVertLevels; ++k)
      d.ru(k, iEdge) = 0.5 * s.u(k, iEdge) * (rho_zz(k, cell1) + rho_zz(k, cell2));
  }
  for (int iCell = 1; iCell <= m.nCells; ++iCell) {  // :7572-7582
    d.rw(1, iCell) = 0.0;
    d.rw(nVertLevels + 1, iCell) = 0.0;
    for (int k = 2; k <= nVertLevels; ++k) {
      d.rw(k, iCell) = s.w(k, iCell) * (m.fzp(k) * rho_zz(k - 1, iCell) + m.fzm(k) * rho_zz(k, iCell)) *
                       (m.fzp(k) * zz(k - 1, iCell) + m.fzm(k) * zz(k, iCell));
    }
  }
  for (int iCell = 1; iCell <= m.nCells; ++iCell) {  // :7587-7600
    for (int i = 1; i <= m.nEdgesOnCell(iCell); ++i) {
      const int iEdge = m.edgesOnCell(i, iCell);
      for (int k = 2; k <= nVertLevels; ++k) {
        const real flux = (m.fzm(k) * d.ru(k, iEdge) + m.fzp(k) * d.ru(k - 1, iEdge));
        d.rw(k, iCell) = d.rw(k, iCell) - m.edgesOnCell_sign(i, iCell) *
                                              (m.zb_cell(k, i, iCell) + std::copysign(1.0, flux) * m.zb3_cell(k, i, iCell)) *
                                              flux * (m.fzp(k) * zz(k - 1, iCell) + m.fzm(k) * zz(k, iCell));
      }
    }
  }
  for (int iCell = 1; iCell <= m.nCells; ++iCell)  // :7605-7609
    for (int k = 1; k <= nVertLevels; ++k) d.rho_p(k, iCell) = rho_zz(k, iCell) - s.rho_base(k, iCell);
  for (int iCell = 1; iCell <= m.nCells; ++iCell)  // :7614-7618
    for (int k = 1; k <= nVertLevels; ++k) d.rtheta_base(k, iCell) = s.theta_base(k, iCell) * s.rho_base(k, iCell);
  for (int iCell = 1; iCell <= m.nCells; ++iCell)  // :7623-7628
    for (int k = 1; k <= nVertLevels; ++k)
      d.rtheta_p(k, iCell) = theta_m(k, iCell) * d.rho_p(k, iCell) +
                             s.rho_base(k, iCell) * (theta_m(k, iCell) - s.theta_base(k, iCell));
  for (int iCell = 1; iCell <= m.nCells; ++iCell) {  // :7633-7638
    auto lane = [&](int k, bool vector) {
      const real xe = zz(k, iCell) * rgas_p0 * (d.rtheta_p(k, iCell) + d.rtheta_base(k, iCell));
      const real xb = zz(k, iCell) * rgas_p0 * (d.rtheta_base(k, iCell));
      d.exner(k, iCell) = vector ? libm::mvec_pow(xe, rcv) : libm::pow(xe, rcv);
      d.exner_base(k, iCell) = vector ? libm::mvec_pow(xb, rcv) : libm::pow(xb, rcv);
    };
    if (nVertLevels == 1) {
      lane(1, false);
      continue;
    }
    int k = 1;
    for (; k + 1 <= nVertLevels; k += 2) {
      lane(k, true);
      lane(k + 1, true);
    }
    if (k <= nVertLevels) lane(k, false);
  }
  for (int iCell = 1; iCell <= m.nCells; ++iCell) {  // :7643-7651
    for (int k = 1; k <= nVertLevels; ++k) {
      d.pressure_p(k, iCell) = zz(k, iCell) * rgas *
                               (d.exner(k, iCell) * d.rtheta_p(k, iCell) +
                                d.rtheta_base(k, iCell) * (d.exner(k, iCell) - d.exner_base(k, iCell)));
      d.pressure_base(k, iCell) = zz(k, iCell) * rgas * d.exner_base(k, iCell) * d.rtheta_base(k, iCell);
    }
  }
}

// mpas_atm_core.F:366-617, in MPAS's order
void atm_mpas_init_block(Mesh& m, State& s, Diag& d, const InitConfig& cfg, std::vector<std::string>& log) {
  if (cfg.config_do_restart) throw Error("atm_mpas_init_block: config_do_restart = true is not ported");
  allocate_init_fields(m, s, d);
  atm_compute_signs(m, log);                                                       // :427
  atm_compute_inverses(m);                                                         // :459-473
  atm_initialize_deformation_weights(m, m.nCells, m.on_a_sphere, m.sphere_radius);  // :478-480
  atm_adv_coef_compression(m);                                                     // :482
  atm_couple_coef_3rd_order(m, cfg);                                               // :484
  // :492-499 ke_vertex / ke_edge are module work arrays of atm_compute_solve_diagnostics (later unit)
  atm_init_coupled_diagnostics(m, s, d);  // :513-521 (cold start)
  // :523-541 atm_compute_solve_diagnostics: later unit
  mpas_initialize_vectors(m);  // :543 mpas_rbf_interp_initialize
  mpas_init_reconstruct(m);    // :544
  // :546-596 first mpas_reconstruct and (moist physics) physics initialisation: later units
  atm_compute_mesh_scaling(m, cfg);   // :599
  atm_compute_damping_coefs(m, cfg);  // :601
  mpas_atm_setup_bdy_masks(m);        // :606
}

}  // namespace cppmpas
