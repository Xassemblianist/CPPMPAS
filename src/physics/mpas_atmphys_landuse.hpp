// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// MPAS-Atmosphere's land-use initialisation and the periodic surface update, statement by statement:
//  * physics/mpas_atmphys_landuse.F (LU): landuse_init_forMPAS (LU:61-321), reading LANDUSE.TBL as gfortran does
//    (the A35 lutype record, list-directed category records; physics/fortran_list_read.hpp);
//  * physics/mpas_atmphys_update_surface.F (US): physics_update_surface (US:52-104);
//  * physics/mpas_atmphys_date_time.F (DT): get_julgmt (DT:51-87), split_date_char (DT:131-149),
//    monthly_interp_to_date (DT:152-225), monthly_min_max (DT:228-258).
// One rank: nCellsSolve = nCells. The arrays are the Registry pools' (io/registry_pools.hpp).
#include <string>

#include "core/mesh.hpp"
#include "io/registry_pools.hpp"

namespace cppmpas::phys {

struct LanduseConfig {
  bool config_do_restart = false, config_frac_seaice = true, config_sfc_albedo = true;
  real xice_threshold = 0;  // mpas_atmphys_manager.F:740-744
};

// LU:61-321. julday: the day of the year of the clock's current time (mpas_atmphys_init.F:181-182). Returns isurban
// (the module variable of mpas_atmphys_landuse); log: MPAS's messages.
int landuse_init_forMPAS(const std::string& table_dir, int julday, const LanduseConfig& cfg, const Mesh& mesh,
                         RegistryPools& p, std::vector<std::string>& log);

// US:52-104. current_date: the clock's date-time string (YYYY-MM-DD_hh:mm:ss).
void physics_update_surface(const std::string& current_date, bool config_sfc_albedo, int nCellsSolve,
                            RegistryPools& p);

// DT:51-87
void get_julgmt(const std::string& date_str, int& julyr, int& julday, real& gmt);

}  // namespace cppmpas::phys
