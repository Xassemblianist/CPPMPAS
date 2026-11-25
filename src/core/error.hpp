// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include <stdexcept>
#include <string>

namespace cppmpas {
// Everything CPPMPAS does not port is an error with a message that names the file, variable or option
// (CONSTITUTION.md 2.4); there is no silent fallback.
struct Error : std::runtime_error {
  explicit Error(const std::string& what) : std::runtime_error(what) {}
};
}  // namespace cppmpas
