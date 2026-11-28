// SPDX-License-Identifier: GPL-3.0-or-later
// See mpas_atmphys_todynamics.hpp. TD = src/core_atmosphere/physics/mpas_atmphys_todynamics.F (MPAS v8.4.2).
#include "physics/mpas_atmphys_todynamics.hpp"

#include <cstddef>

#include "core/backend.hpp"
#include "physics/mpas_atmphys_interface.hpp"

#include "core/error.hpp"
#include "physics/mpas_atmphys_constants.hpp"

namespace cppmpas::phys {

namespace {

// TD:501-552 over every edge 1..nEdges (one thread per edge)
void tend_toEdges(const GetTendArrays& m, FView<const real> Ux_tend, FView<const real> Uy_tend, FView<real> U_tend) {
  parallel_for(1, m.nEdges, [=] CPPMPAS_LAMBDA(int iEdge) {
    const auto& n = m.edgeNormalVectors;
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
  });
}

}  // namespace

// One kernel per loop of TD (one thread per cell or edge, the k and scalar loops inside the thread).
void physics_get_tend(GetTendPools& mp, int rk_step, int dynamics_substep, FView<real> tend_ru_physics,
                      FView<real> tend_rtheta_physics, FView<real> tend_rho_physics,
                      const std::function<void(const char*)>& exchange_halo_group) {
  using namespace physics_constants;
  const std::string& pbl = mp.pbl_scheme;
  const std::string& cu = mp.convection_scheme;
  if (pbl != "off" && pbl != "bl_ysu" && pbl != "bl_mynn")
    throw Error("physics_get_tend: config_pbl_scheme = '" + pbl + "' has no verified port");
  if (pbl == "bl_mynn" && mp.microp_scheme == "mp_thompson_aerosols")  // TD:414-422
    throw Error("physics_get_tend: bl_mynn with mp_thompson_aerosols (nc / nifa / nwfa tendencies) has no verified port");
  if (pbl == "bl_mynn" && (mp.index_qs <= 0 || mp.index_ni <= 0))  // -1: an inactive constituent
    throw Error("physics_get_tend: bl_mynn needs the scalars qs and ni (TD:409-410)");
  if (pbl == "bl_mynn" && (!mp.rqsblten.p || !mp.rniblten.p))
    throw Error("physics_get_tend: bl_mynn's rqsblten / rniblten (tend_physics, bl_mynn_in) are not set");
  const GetTendArrays m = mp;
  const int nV = m.nVertLevels, nCells = m.nCells, nEdges = m.nEdges;
  const int nCellsSolve = m.nCellsSolve, nEdgesSolve = m.nEdgesSolve;
  const int num_scalars = m.tend_scalars.e0;
  const auto tend_scalars = m.tend_scalars;
  const auto mass = m.mass;
  const auto mass_edge = m.mass_edge;

  // TD:210-216
  WorkArray tend_th_a(static_cast<std::size_t>(nV) * nCellsSolve, 0.);
  real* tend_th_v = tend_th_a.data();
#define tend_th(k, i) tend_th_v[static_cast<std::size_t>((i) - 1) * nV + ((k) - 1)]
  parallel_for(1, nCells + 1, [=] CPPMPAS_LAMBDA(int i) {
    for (int k = 1; k <= nV; ++k)
      for (int s = 1; s <= num_scalars; ++s) tend_scalars(s, k, i) = 0.;
  });
  parallel_for(1, nEdges + 1, [=] CPPMPAS_LAMBDA(int i) {
    for (int k = 1; k <= nV; ++k) tend_ru_physics(k, i) = 0.;
  });
  parallel_for(1, nCells + 1, [=] CPPMPAS_LAMBDA(int i) {
    for (int k = 1; k <= nV; ++k) {
      tend_rtheta_physics(k, i) = 0.;
      tend_rho_physics(k, i) = 0.;
    }
  });
  // physics_get_tend_work (TD:304-498); tend_u = tend_ru_physics, tend_theta = tend_rtheta_physics
  const auto tend_u = tend_ru_physics;
  const auto tend_theta = tend_rtheta_physics;

  if (pbl != "off") {  // TD:382-426
    if (rk_step == 1 && dynamics_substep == 1) {
      exchange_halo_group("physics:blten");
      tend_toEdges(m, m.rublten, m.rvblten, m.rublten_Edge);
      parallel_for(1, nEdges, [=] CPPMPAS_LAMBDA(int i) {
        for (int k = 1; k <= nV; ++k) m.tend_u_phys(k, i) = m.rublten_Edge(k, i);
      });
    }
    parallel_for(1, nEdgesSolve, [=] CPPMPAS_LAMBDA(int i) {
      for (int k = 1; k <= nV; ++k) tend_u(k, i) = tend_u(k, i) + m.rublten_Edge(k, i) * mass_edge(k, i);
    });
    parallel_for(1, nCellsSolve, [=] CPPMPAS_LAMBDA(int i) {
      for (int k = 1; k <= nV; ++k) {
        tend_th(k, i) = tend_th(k, i) + m.rthblten(k, i) * mass(k, i);
        tend_scalars(m.index_qv, k, i) = tend_scalars(m.index_qv, k, i) + m.rqvblten(k, i) * mass(k, i);
        tend_scalars(m.index_qc, k, i) = tend_scalars(m.index_qc, k, i) + m.rqcblten(k, i) * mass(k, i);
        tend_scalars(m.index_qi, k, i) = tend_scalars(m.index_qi, k, i) + m.rqiblten(k, i) * mass(k, i);
      }
    });
    if (pbl == "bl_mynn")  // TD:405-412
      parallel_for(1, nCellsSolve, [=] CPPMPAS_LAMBDA(int i) {
        for (int k = 1; k <= nV; ++k) {
          tend_scalars(m.index_qs, k, i) = tend_scalars(m.index_qs, k, i) + m.rqsblten(k, i) * mass(k, i);
          tend_scalars(m.index_ni, k, i) = tend_scalars(m.index_ni, k, i) + m.rniblten(k, i) * mass(k, i);
        }
      });
  }

  if (cu != "off") {  // TD:430-465
    parallel_for(1, nCellsSolve, [=] CPPMPAS_LAMBDA(int i) {
      for (int k = 1; k <= nV; ++k) {
        tend_th(k, i) = tend_th(k, i) + m.rthcuten(k, i) * mass(k, i);
        tend_scalars(m.index_qv, k, i) = tend_scalars(m.index_qv, k, i) + m.rqvcuten(k, i) * mass(k, i);
        tend_scalars(m.index_qc, k, i) = tend_scalars(m.index_qc, k, i) + m.rqccuten(k, i) * mass(k, i);
        tend_scalars(m.index_qi, k, i) = tend_scalars(m.index_qi, k, i) + m.rqicuten(k, i) * mass(k, i);
      }
    });
    if (cu == "cu_kain_fritsch") {  // TD:441-447
      parallel_for(1, nCellsSolve, [=] CPPMPAS_LAMBDA(int i) {
        for (int k = 1; k <= nV; ++k) {
          tend_scalars(m.index_qr, k, i) = tend_scalars(m.index_qr, k, i) + m.rqrcuten(k, i) * mass(k, i);
          tend_scalars(m.index_qs, k, i) = tend_scalars(m.index_qs, k, i) + m.rqscuten(k, i) * mass(k, i);
        }
      });
    } else if (cu == "cu_tiedtke" || cu == "cu_ntiedtke") {  // TD:449-461
      if (rk_step == 1 && dynamics_substep == 1) {
        exchange_halo_group("physics:cuten");
        tend_toEdges(m, m.rucuten, m.rvcuten, m.rucuten_Edge);
        parallel_for(1, nEdges, [=] CPPMPAS_LAMBDA(int i) {
          for (int k = 1; k <= nV; ++k) m.tend_u_phys(k, i) = m.tend_u_phys(k, i) + m.rucuten_Edge(k, i);
        });
      }
      parallel_for(1, nEdgesSolve, [=] CPPMPAS_LAMBDA(int i) {
        for (int k = 1; k <= nV; ++k) tend_u(k, i) = tend_u(k, i) + m.rucuten_Edge(k, i) * mass_edge(k, i);
      });
    }
  }

  if (mp.radt_lw_scheme != "off")  // TD:469-475
    parallel_for(1, nCellsSolve, [=] CPPMPAS_LAMBDA(int i) {
      for (int k = 1; k <= nV; ++k) tend_th(k, i) = tend_th(k, i) + m.rthratenlw(k, i) * mass(k, i);
    });
  if (mp.radt_sw_scheme != "off")  // TD:479-485
    parallel_for(1, nCellsSolve, [=] CPPMPAS_LAMBDA(int i) {
      for (int k = 1; k <= nV; ++k) tend_th(k, i) = tend_th(k, i) + m.rthratensw(k, i) * mass(k, i);
    });

  parallel_for(1, nCellsSolve, [=] CPPMPAS_LAMBDA(int i) {  // TD:489-495
    for (int k = 1; k <= nV; ++k) {
      const real coeff = (1. + R_v / R_d * m.scalars(m.index_qv, k, i));
      tend_th(k, i) = coeff * tend_th(k, i) + R_v / R_d * m.theta_m(k, i) * tend_scalars(m.index_qv, k, i) / coeff;
      tend_theta(k, i) = tend_theta(k, i) + tend_th(k, i);
    }
  });
#undef tend_th
}

}  // namespace cppmpas::phys
