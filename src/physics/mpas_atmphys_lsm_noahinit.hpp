// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// mpas_atmphys_lsm_noahinit.F (MPAS v8.4.2): soil_veg_gen_parm (:210-589), the reads of VEGPARM.TBL, SOILPARM.TBL and
// GENPARM.TBL into module_sf_noahlsm's tables, with gfortran's list-directed input (one rank: the DM_BCAST calls
// copy nothing). MPAS's own reader, not WRF's SOIL_VEG_GEN_PARM (module_sf_noahdrv.F:1991-2332; CPPWRF's
// noah_tables.cpp reads the tables its own way): MPAS skips the three records before LCZ_1 where WRF reads one as
// text to refuse old tables, and MPAS's SOILPARM.TBL size check is an empty IF where WRF stops.
#include <string>

#include "physics/sf_noah.hpp"

namespace cppmpas::phys {

// soil_veg_gen_parm(dminfo, mminlu, mminsl): table_dir holds the three tables (MPAS opens them in the run
// directory). Every FATAL_ERROR of the routine, and every read the Fortran would stop on, is an Error.
void soil_veg_gen_parm(const std::string& table_dir, const std::string& mminlu, const std::string& mminsl,
                       noah::NoahTables& t);

}  // namespace cppmpas::phys
