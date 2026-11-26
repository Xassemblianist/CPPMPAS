// SPDX-License-Identifier: GPL-3.0-or-later
// See mpas_atmphys_todynamics.hpp. TD = src/core_atmosphere/physics/mpas_atmphys_todynamics.F (MPAS v8.4.2).
#include "physics/mpas_atmphys_todynamics.hpp"

#include <vector>

#include "core/error.hpp"
#include "physics/mpas_atmphys_constants.hpp"

namespace cppmpas::phys {

namespace {

// TD:501-552 over every edge 1..nEdges
void tend_toEdges(const GetTendPools& m, FView<const real> Ux_tend, FView<const real> Uy_tend, FView<real> U_tend) {
  const auto& n = m.edgeNormalVectors;
  for (int iEdge = 1; iEdge <= m.nEdges; ++iEdge) {
    const int cell1 = m.cellsOnEdge(1, iEdge);
    const int cell2 = m.cellsOnEdge(2, iEdge);
    for (int k = 1; k <= m.nVertLevels; ++k)  // TD:538-549
      U_tend(k, iEdge) = Ux_tend(k, cell1) * 0.5 * (n(1, iEdge) * m.east(1, cell1)
                                                  + n(2, iEdge) * m.east(2, cell1)
                                                  + n(3, iEdge) * m.east(3, cell1))
                       + Uy_tend(k, cell1) * 0.5 * (n(1, iEdge) * m.north(1, cell1)
                                                  + n(2, iEdge) * m.north(2, cell1)
                                                  + n(3, iEdge) * m.north(3, cell1))
                       + Ux_tend(k, cell2) * 0.5 * (n(1, iEdge) * m.east(1, cell2)
                                                  + n(2, iEdge) * m.east(2, cell2)
                                                  + n(3, iEdge) * m.east(3, cell2))
                       + Uy_tend(k, cell2) * 0.5 * (n(1, iEdge) * m.north(1, cell2)
                                                  + n(2, iEdge) * m.north(2, cell2)
                                                  + n(3, iEdge) * m.north(3, cell2));
  }
}

}  // namespace

void physics_get_tend(GetTendPools& m, int rk_step, int dynamics_substep, FView<real> tend_ru_physics,
                      FView<real> tend_rtheta_physics, FView<real> tend_rho_physics,
                      const std::function<void(const char*)>& exchange_halo_group) {
  using namespace physics_constants;
  const std::string& pbl = m.pbl_scheme;
  const std::string& cu = m.convection_scheme;
  if (pbl != "off" && pbl != "bl_ysu")  // TD:405-425: bl_mynn's qs / ni / nc / nifa / nwfa tendencies
    throw Error("physics_get_tend: config_pbl_scheme = '" + pbl + "' has no verified port");
  const int nV = m.nVertLevels, nCells = m.nCells, nEdges = m.nEdges;
  const int nCellsSolve = m.nCellsSolve, nEdgesSolve = m.nEdgesSolve;
  const int num_scalars = m.tend_scalars.e0;
  auto& tend_scalars = m.tend_scalars;
  const auto& mass = m.mass;
  const auto& mass_edge = m.mass_edge;

  // TD:210-216
  std::vector<real> tend_th_v(static_cast<std::size_t>(nV) * nCellsSolve, 0.);
  auto tend_th = [&](int k, int i) -> real& { return tend_th_v[static_cast<std::size_t>(i - 1) * nV + (k - 1)]; };
  for (int i = 1; i <= nCells + 1; ++i)
    for (int k = 1; k <= nV; ++k)
      for (int s = 1; s <= num_scalars; ++s) tend_scalars(s, k, i) = 0.;
  for (int i = 1; i <= nEdges + 1; ++i)
    for (int k = 1; k <= nV; ++k) tend_ru_physics(k, i) = 0.;
  for (int i = 1; i <= nCells + 1; ++i)
    for (int k = 1; k <= nV; ++k) {
      tend_rtheta_physics(k, i) = 0.;
      tend_rho_physics(k, i) = 0.;
    }
  // physics_get_tend_work (TD:304-498); tend_u = tend_ru_physics, tend_theta = tend_rtheta_physics
  auto& tend_u = tend_ru_physics;
  auto& tend_theta = tend_rtheta_physics;

  if (pbl != "off") {  // TD:382-426
    if (rk_step == 1 && dynamics_substep == 1) {
      exchange_halo_group("physics:blten");
      tend_toEdges(m, m.rublten, m.rvblten, m.rublten_Edge);
      for (int i = 1; i <= nEdges; ++i)
        for (int k = 1; k <= nV; ++k) m.tend_u_phys(k, i) = m.rublten_Edge(k, i);
    }
    for (int i = 1; i <= nEdgesSolve; ++i)
      for (int k = 1; k <= nV; ++k) tend_u(k, i) = tend_u(k, i) + m.rublten_Edge(k, i) * mass_edge(k, i);
    for (int i = 1; i <= nCellsSolve; ++i)
      for (int k = 1; k <= nV; ++k) {
        tend_th(k, i) = tend_th(k, i) + m.rthblten(k, i) * mass(k, i);
        tend_scalars(m.index_qv, k, i) = tend_scalars(m.index_qv, k, i) + m.rqvblten(k, i) * mass(k, i);
        tend_scalars(m.index_qc, k, i) = tend_scalars(m.index_qc, k, i) + m.rqcblten(k, i) * mass(k, i);
        tend_scalars(m.index_qi, k, i) = tend_scalars(m.index_qi, k, i) + m.rqiblten(k, i) * mass(k, i);
      }
  }

  if (cu != "off") {  // TD:430-465
    for (int i = 1; i <= nCellsSolve; ++i)
      for (int k = 1; k <= nV; ++k) {
        tend_th(k, i) = tend_th(k, i) + m.rthcuten(k, i) * mass(k, i);
        tend_scalars(m.index_qv, k, i) = tend_scalars(m.index_qv, k, i) + m.rqvcuten(k, i) * mass(k, i);
        tend_scalars(m.index_qc, k, i) = tend_scalars(m.index_qc, k, i) + m.rqccuten(k, i) * mass(k, i);
        tend_scalars(m.index_qi, k, i) = tend_scalars(m.index_qi, k, i) + m.rqicuten(k, i) * mass(k, i);
      }
    if (cu == "cu_kain_fritsch") {  // TD:441-447
      for (int i = 1; i <= nCellsSolve; ++i)
        for (int k = 1; k <= nV; ++k) {
          tend_scalars(m.index_qr, k, i) = tend_scalars(m.index_qr, k, i) + m.rqrcuten(k, i) * mass(k, i);
          tend_scalars(m.index_qs, k, i) = tend_scalars(m.index_qs, k, i) + m.rqscuten(k, i) * mass(k, i);
        }
    } else if (cu == "cu_tiedtke" || cu == "cu_ntiedtke") {  // TD:449-461
      if (rk_step == 1 && dynamics_substep == 1) {
        exchange_halo_group("physics:cuten");
        tend_toEdges(m, m.rucuten, m.rvcuten, m.rucuten_Edge);
        for (int i = 1; i <= nEdges; ++i)
          for (int k = 1; k <= nV; ++k) m.tend_u_phys(k, i) = m.tend_u_phys(k, i) + m.rucuten_Edge(k, i);
      }
      for (int i = 1; i <= nEdgesSolve; ++i)
        for (int k = 1; k <= nV; ++k) tend_u(k, i) = tend_u(k, i) + m.rucuten_Edge(k, i) * mass_edge(k, i);
    }
  }

  if (m.radt_lw_scheme != "off")  // TD:469-475
    for (int i = 1; i <= nCellsSolve; ++i)
      for (int k = 1; k <= nV; ++k) tend_th(k, i) = tend_th(k, i) + m.rthratenlw(k, i) * mass(k, i);
  if (m.radt_sw_scheme != "off")  // TD:479-485
    for (int i = 1; i <= nCellsSolve; ++i)
      for (int k = 1; k <= nV; ++k) tend_th(k, i) = tend_th(k, i) + m.rthratensw(k, i) * mass(k, i);

  for (int i = 1; i <= nCellsSolve; ++i)  // TD:489-495
    for (int k = 1; k <= nV; ++k) {
      const real coeff = (1. + R_v / R_d * m.scalars(m.index_qv, k, i));
      tend_th(k, i) = coeff * tend_th(k, i) + R_v / R_d * m.theta_m(k, i) * tend_scalars(m.index_qv, k, i) / coeff;
      tend_theta(k, i) = tend_theta(k, i) + tend_th(k, i);
    }
}

}  // namespace cppmpas::phys
