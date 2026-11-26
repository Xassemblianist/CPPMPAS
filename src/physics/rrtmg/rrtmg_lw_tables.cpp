// SPDX-License-Identifier: GPL-3.0-or-later AND LicenseRef-AER-RRTMG
// Shared with CPPWRF (SHARED_WITH_CPPWRF.md): CPPWRF 7fac5bf src/phys/ra_rrtmg_lw/rrtmg_lw_tables.cpp, copied unchanged except the namespace
// (cppmpas::phys::rrtmg), the CPPMPAS_ macro prefix, the include paths and the real type.
#include "physics/rrtmg/rrtmg_lw_tables.hpp"

#include <cstddef>
#include <iostream>

#include "physics/rrtmg/table_dump.hpp"

namespace cppmpas::phys::rrtmg {

namespace {

// the int32 members of RrtmgLwTablesOf<R>: ngc .. ngm, one run before wt
template <class R>
constexpr std::array<TableBytes, 1> lw_table_ints() {
  using T = RrtmgLwTablesOf<R>;
  return {TableBytes{offsetof(T, ngc), offsetof(T, wt)}};
}

}  // namespace

bool load_rrtmg_lw_tables(const std::string& filepath, RrtmgLwTables& tables) {
  if (read_table_dump<RrtmgLwTablesOf<double>>(filepath, tables, lw_table_ints<double>(), lw_table_ints<rb>()))
    return true;
  std::cerr << "Error: could not read rrtmg_lw_tables.bin from " << filepath << "\n";
  return false;
}

}  // namespace cppmpas::phys::rrtmg
