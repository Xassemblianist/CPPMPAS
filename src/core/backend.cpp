// SPDX-License-Identifier: GPL-3.0-or-later
// Host side of core/backend.hpp shared by the CPU and the CUDA build.
#include "core/backend.hpp"

namespace cppmpas {

namespace detail {
void throw_kernel_error(const char* head, bool has_value, long value, const char* tail) {
  std::string msg = head ? head : "";
  if (has_value) msg += std::to_string(value);
  if (tail) msg += tail;
  throw Error(msg);
}
}  // namespace detail

void fits_column(int nVertLevels, const char* routine) {
#if !defined(CPPMPAS_USE_CUDA)
  (void)nVertLevels;
  (void)routine;
#else
  if (nVertLevels + 2 > kMaxColumn)
    throw Error(std::string(routine) + ": nVertLevels = " + std::to_string(nVertLevels) + " is more than the " +
                std::to_string(kMaxColumn - 2) + " levels of the GPU's column work arrays (core/backend.hpp kMaxColumn)");
#endif
}

}  // namespace cppmpas
