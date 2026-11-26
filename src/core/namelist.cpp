// SPDX-License-Identifier: GPL-3.0-or-later
#include "core/namelist.hpp"

#include <cctype>
#include <cstdlib>
#include <fstream>
#include <sstream>

#include "core/error.hpp"

namespace cppmpas {

namespace {

std::string lower(std::string s) {
  for (auto& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  return s;
}
std::string trim(const std::string& s) {
  std::size_t a = 0, b = s.size();
  while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
  while (b > a && (std::isspace(static_cast<unsigned char>(s[b - 1])) || s[b - 1] == ',')) --b;
  return s.substr(a, b - a);
}

// a real as gfortran's list-directed namelist read converts it: correctly rounded (strtod); d/D exponents
real parse_real(const std::string& key, std::string v) {
  for (auto& c : v)
    if (c == 'd' || c == 'D') c = 'e';
  char* end = nullptr;
  const double x = std::strtod(v.c_str(), &end);
  if (v.empty() || *end != '\0') throw Error("namelist: " + key + " = '" + v + "' is not a real number");
  return x;
}
bool parse_logical(const std::string& key, const std::string& v) {
  std::string s = lower(v);
  if (!s.empty() && s[0] == '.') s = s.substr(1);
  if (!s.empty() && (s[0] == 't')) return true;
  if (!s.empty() && (s[0] == 'f')) return false;
  throw Error("namelist: " + key + " = '" + v + "' is not a logical");
}

}  // namespace

std::map<std::string, std::string> read_namelist_file(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw Error("cannot open namelist file " + path);
  std::map<std::string, std::string> kv;
  std::string line;
  bool in_group = false;
  while (std::getline(in, line)) {
    // strip a comment that is not inside quotes
    bool q1 = false, q2 = false;
    for (std::size_t i = 0; i < line.size(); ++i) {
      if (line[i] == '\'' && !q2) q1 = !q1;
      else if (line[i] == '"' && !q1) q2 = !q2;
      else if (line[i] == '!' && !q1 && !q2) { line.resize(i); break; }
    }
    const std::string t = trim(line);
    if (t.empty()) continue;
    if (t[0] == '&') { in_group = true; continue; }
    if (t == "/") { in_group = false; continue; }
    if (!in_group) throw Error(path + ": text outside a namelist group: " + t);
    const auto eq = t.find('=');
    if (eq == std::string::npos) throw Error(path + ": not a 'key = value' line: " + t);
    const std::string key = lower(trim(t.substr(0, eq)));
    std::string val = trim(t.substr(eq + 1));
    if (val.size() >= 2 && (val.front() == '\'' || val.front() == '"') && val.back() == val.front())
      val = val.substr(1, val.size() - 2);
    if (kv.count(key)) throw Error(path + ": option " + key + " is set twice");
    kv[key] = val;
  }
  return kv;
}

InitConfig read_init_config(const std::string& path) {
  InitConfig c;
  if (path.empty()) return c;
  const auto kv = read_namelist_file(path);
  auto r = [&](const char* k, real& x) {
    const auto it = kv.find(k);
    if (it != kv.end()) x = parse_real(k, it->second);
  };
  auto l = [&](const char* k, bool& x) {
    const auto it = kv.find(k);
    if (it != kv.end()) x = parse_logical(k, it->second);
  };
  r("config_dt", c.config_dt);
  r("config_apvm_upwinding", c.config_apvm_upwinding);
  r("config_coef_3rd_order", c.config_coef_3rd_order);
  l("config_h_scalewithmesh", c.config_h_ScaleWithMesh);
  r("config_zd", c.config_zd);
  r("config_xnutr", c.config_xnutr);
  r("config_epssm_minimum", c.config_epssm_minimum);
  r("config_epssm_maximum", c.config_epssm_maximum);
  r("config_epssm_transition_bottom_z", c.config_epssm_transition_bottom_z);
  r("config_epssm_transition_top_z", c.config_epssm_transition_top_z);
  l("config_do_restart", c.config_do_restart);
  l("config_do_dacycling", c.config_do_DAcycling);
  if (c.config_do_restart)
    throw Error(path + ": config_do_restart = true: the restart start-up (mpas_atm_core.F:507-523, fields read from "
                "the restart stream) is not ported");
  return c;
}

}  // namespace cppmpas
