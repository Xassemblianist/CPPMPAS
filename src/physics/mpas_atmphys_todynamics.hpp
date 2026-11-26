// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// MPAS-Atmosphere's physics tendencies for the dynamics (physics/mpas_atmphys_todynamics.F, TD): physics_get_tend
// (TD:98-279) with physics_get_tend_work (TD:304-498) and tend_toEdges (TD:501-552), statement by statement.
//
// Ported: config_pbl_scheme 'off' or 'bl_ysu' (TD:382-403; 'bl_mynn', TD:405-425, is refused), every
// config_convection_scheme branch (TD:430-465: the common th / qv / qc / qi terms, cu_kain_fritsch's qr / qs,
// cu_tiedtke / cu_ntiedtke's momentum through tend_toEdges), the longwave and shortwave heating (TD:469-485) and the
// moist theta coupling (TD:489-495). The halo exchanges ('physics:blten', 'physics:cuten') are the caller's
// exchange_halo_group.
//
// Reference object code (bin/atmosphere_model, sha256 524a1cf4...): physics_get_tend (physics_get_tend_work inlined)
// and tend_toedges.constprop.0 have no FMA, no ymm register and no libm call. MPAS's file compiled standalone with the
// reference flags has the reference's arithmetic instruction counts (tend_toedges: 379 instructions, the set-up of 14
// differs: the stub pool calls).
#include <functional>
#include <string>

#include "core/field.hpp"
#include "core/precision.hpp"

namespace cppmpas::phys {

struct GetTendPools {
  int nCells = 0, nEdges = 0, nCellsSolve = 0, nEdgesSolve = 0, nVertLevels = 0;
  // configs
  std::string pbl_scheme, convection_scheme, microp_scheme, radt_lw_scheme, radt_sw_scheme;
  // mesh (tend_toEdges)
  FView<const int> cellsOnEdge;
  FView<const real> edgeNormalVectors, east, north;
  // state: theta_m and scalars of time level 1, rho_zz of time level 2 (TD:163-165); diag rho_edge (TD:166)
  FView<const real> theta_m, scalars, mass, mass_edge;
  int index_qv = 0, index_qc = 0, index_qr = 0, index_qi = 0, index_qs = 0;
  // diag
  FView<real> tend_u_phys;
  // tend_physics
  FView<const real> rublten, rvblten, rthblten, rqvblten, rqcblten, rqiblten;
  FView<const real> rucuten, rvcuten, rthcuten, rqvcuten, rqccuten, rqrcuten, rqicuten, rqscuten;
  FView<const real> rthratenlw, rthratensw;
  FView<real> rublten_Edge, rucuten_Edge;
  // tend
  FView<real> tend_scalars;
};

// TD:98-279. tend_ru_physics (nVertLevels, nEdges+1), tend_rtheta_physics and tend_rho_physics (nVertLevels, nCells+1)
// are zeroed whole and filled, as the dynamics' arrays of mpas_atm_time_integration.F.
void physics_get_tend(GetTendPools& m, int rk_step, int dynamics_substep, FView<real> tend_ru_physics,
                      FView<real> tend_rtheta_physics, FView<real> tend_rho_physics,
                      const std::function<void(const char*)>& exchange_halo_group);

}  // namespace cppmpas::phys
