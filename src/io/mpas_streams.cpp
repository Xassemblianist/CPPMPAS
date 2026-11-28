// SPDX-License-Identifier: GPL-3.0-or-later
// MPAS's stream manager for the output streams of a run: see mpas_streams.hpp.
#include "io/mpas_streams.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <map>
#include <set>
#include <sstream>

#include "core/error.hpp"
#include "io/registry.hpp"

namespace cppmpas {

namespace {

std::string slurp(const std::string& path) {
  std::ifstream f(path);
  if (!f) throw Error("streams: cannot read " + path);
  std::stringstream ss;
  ss << f.rdbuf();
  return ss.str();
}

struct Tag {
  std::string name;
  bool closing = false, self_closing = false;
  std::map<std::string, std::string> attrs;
};

// The tags of the XML file, comments skipped (xml_stream_parser.c reads it with ezxml).
std::vector<Tag> tags_of(const std::string& text, const std::string& path) {
  std::vector<Tag> out;
  std::size_t i = 0;
  while ((i = text.find('<', i)) != std::string::npos) {
    if (text.compare(i, 4, "<!--") == 0) {
      const std::size_t e = text.find("-->", i);
      if (e == std::string::npos) throw Error("streams: unterminated comment in " + path);
      i = e + 3;
      continue;
    }
    if (text.compare(i, 2, "<?") == 0) {
      i = text.find("?>", i);
      if (i == std::string::npos) throw Error("streams: bad declaration in " + path);
      continue;
    }
    const std::size_t e = text.find('>', i);
    if (e == std::string::npos) throw Error("streams: unterminated tag in " + path);
    std::string body = text.substr(i + 1, e - i - 1);
    i = e + 1;
    Tag t;
    if (!body.empty() && body[0] == '/') {
      t.closing = true;
      body.erase(0, 1);
    }
    if (!body.empty() && body.back() == '/') {
      t.self_closing = true;
      body.pop_back();
    }
    std::size_t p = 0;
    while (p < body.size() && !std::isspace(static_cast<unsigned char>(body[p]))) ++p;
    t.name = body.substr(0, p);
    while (p < body.size()) {
      while (p < body.size() && std::isspace(static_cast<unsigned char>(body[p]))) ++p;
      if (p >= body.size()) break;
      const std::size_t eq = body.find('=', p);
      if (eq == std::string::npos) throw Error("streams: bad attribute in <" + t.name + "> of " + path);
      std::string key = body.substr(p, eq - p);
      while (!key.empty() && std::isspace(static_cast<unsigned char>(key.back()))) key.pop_back();
      std::size_t q = eq + 1;
      while (q < body.size() && std::isspace(static_cast<unsigned char>(body[q]))) ++q;
      if (q >= body.size() || (body[q] != '"' && body[q] != '\''))
        throw Error("streams: unquoted attribute " + key + " in " + path);
      const char quote = body[q];
      const std::size_t qe = body.find(quote, q + 1);
      if (qe == std::string::npos) throw Error("streams: unterminated attribute " + key + " in " + path);
      if (!t.attrs.emplace(key, body.substr(q + 1, qe - q - 1)).second)
        throw Error("streams: attribute " + key + " given twice in " + path);
      p = qe + 1;
    }
    out.push_back(std::move(t));
  }
  return out;
}

// the field names of a <file> list: one name per line (xml_stream_parser.c:parse_xml_tag_name ... fgets per line,
// blanks trimmed)
std::vector<std::string> read_field_list(const std::string& path) {
  std::ifstream f(path);
  if (!f) throw Error("streams: cannot read the field list " + path);
  std::vector<std::string> out;
  std::string line;
  while (std::getline(f, line)) {
    const auto b = line.find_first_not_of(" \t\r");
    if (b == std::string::npos) continue;
    const auto e = line.find_last_not_of(" \t\r");
    out.push_back(line.substr(b, e - b + 1));
  }
  return out;
}

bool any_active(const std::string& packages, const ActivePackages& ap) {
  if (packages.empty()) return true;
  std::stringstream ss(packages);
  std::string p;
  while (std::getline(ss, p, ';'))
    if (ap.on.count(p)) return true;
  return false;
}

}  // namespace

std::vector<StreamSpec> read_streams_file(const std::string& path) {
  const std::string dir = path.find('/') == std::string::npos ? "." : path.substr(0, path.rfind('/'));
  const auto tags = tags_of(slurp(path), path);
  std::vector<StreamSpec> out;
  StreamSpec* cur = nullptr;
  bool in_streams = false;
  static const std::set<std::string> kKnown = {"name", "type", "filename_template", "filename_interval",
                                               "input_interval", "output_interval", "packages"};
  for (const auto& t : tags) {
    if (t.name == "streams") {
      in_streams = !t.closing;
      continue;
    }
    if (!in_streams) throw Error("streams: <" + t.name + "> outside <streams> in " + path);
    if (t.name == "immutable_stream" || t.name == "stream") {
      if (t.closing) {
        cur = nullptr;
        continue;
      }
      StreamSpec s;
      s.immutable = t.name == "immutable_stream";
      for (const auto& [k, v] : t.attrs)
        if (!kKnown.count(k))
          throw Error("streams: the stream attribute " + k + " (stream '" +
                      (t.attrs.count("name") ? t.attrs.at("name") : "?") + "') has no verified port");
      auto get = [&](const char* k) { return t.attrs.count(k) ? t.attrs.at(k) : std::string(); };
      s.name = get("name");
      if (s.name.empty()) throw Error("streams: a stream without a name in " + path);
      const std::string type = get("type");
      if (type == "input") s.input = true;
      else if (type == "output") s.output = true;
      else if (type == "input;output" || type == "output;input") s.input = s.output = true;
      else if (type != "none") throw Error("streams: stream '" + s.name + "' has type '" + type + "'");
      s.filename_template = get("filename_template");
      s.filename_interval = get("filename_interval");
      if (t.attrs.count("input_interval")) s.input_interval = get("input_interval");
      if (t.attrs.count("output_interval")) s.output_interval = get("output_interval");
      s.packages = get("packages");
      if (s.immutable) {
        if (s.name == "input") s.fields = input_stream_vars();
        else if (s.name == "restart") s.fields = restart_stream_vars();
        // the other immutable streams (iau, lbc_in, ugwp_*, da_state) hold Registry fields CPPMPAS does not list;
        // they must be inactive or never ring (checked in make_stream_manager)
      }
      out.push_back(std::move(s));
      cur = t.self_closing ? nullptr : &out.back();
      continue;
    }
    if (t.closing) continue;
    if (!cur) throw Error("streams: <" + t.name + "> outside a stream in " + path);
    if (cur->immutable) throw Error("streams: immutable stream '" + cur->name + "' lists fields");
    if (t.name == "var" || t.name == "var_array") {
      if (t.attrs.size() != 1 || !t.attrs.count("name"))
        throw Error("streams: <" + t.name + "> of stream '" + cur->name + "' needs exactly a name");
      cur->fields.push_back(t.attrs.at("name"));
    } else if (t.name == "file") {
      if (t.attrs.size() != 1 || !t.attrs.count("name"))
        throw Error("streams: <file> of stream '" + cur->name + "' needs exactly a name");
      for (auto& f : read_field_list(dir + "/" + t.attrs.at("name"))) cur->fields.push_back(f);
    } else {
      throw Error("streams: <" + t.name + "> in stream '" + cur->name + "' has no verified port");
    }
  }
  return out;
}

long long interval_seconds(const std::string& s) {
  // mpas_timekeeping.F mpas_set_timeInterval(timeString): [D_]h:m:s, [D_]m:s? -- MPAS reads "DD_hh:mm:ss",
  // "hh:mm:ss", "mm:ss" and "ss" (mpas_timekeeping.F:1435-1520)
  auto num = [&](const std::string& x) {
    if (x.empty() || x.find_first_not_of("0123456789") != std::string::npos)
      throw Error("time interval '" + s + "' is not of the form [D_]hh:mm:ss");
    return std::stoll(x);
  };
  if (s.find('.') != std::string::npos) throw Error("time interval '" + s + "': fractions have no verified port");
  long long days = 0;
  std::string rest = s;
  if (const auto u = s.find('_'); u != std::string::npos) {
    days = num(s.substr(0, u));
    rest = s.substr(u + 1);
  }
  std::vector<std::string> parts;
  std::stringstream ss(rest);
  std::string p;
  while (std::getline(ss, p, ':')) parts.push_back(p);
  if (parts.empty() || parts.size() > 3) throw Error("time interval '" + s + "' is not of the form [D_]hh:mm:ss");
  long long h = 0, m = 0, sec = 0;
  if (parts.size() == 3) {
    h = num(parts[0]);
    m = num(parts[1]);
    sec = num(parts[2]);
  } else if (parts.size() == 2) {
    m = num(parts[0]);
    sec = num(parts[1]);
  } else {
    sec = num(parts[0]);
  }
  return ((days * 24 + h) * 60 + m) * 60 + sec;
}

std::string expand_filename(const std::string& tmpl, const std::string& stamp) {
  if (stamp.size() < 19) throw Error("expand_filename: bad time stamp '" + stamp + "'");
  std::string out;
  for (std::size_t i = 0; i < tmpl.size(); ++i) {
    if (tmpl[i] != '$' || i + 1 >= tmpl.size()) {
      out += tmpl[i];
      continue;
    }
    const char c = tmpl[++i];
    switch (c) {
      case 'Y': out += stamp.substr(0, 4); break;
      case 'M': out += stamp.substr(5, 2); break;
      case 'D': out += stamp.substr(8, 2); break;
      case 'h': out += stamp.substr(11, 2); break;
      case 'm': out += stamp.substr(14, 2); break;
      case 's': out += stamp.substr(17, 2); break;
      default: throw Error("filename_template '" + tmpl + "': $" + std::string(1, c) + " has no verified port");
    }
  }
  return out;
}

bool StreamManager::field_will_be_written(const std::string& field, long long now) const {
  for (const auto& s : outputs)
    if (s.ringing(now) && std::find(s.fields.begin(), s.fields.end(), field) != s.fields.end()) return true;
  return false;
}

int StreamManager::output_inclusion_count(const std::string& field) const {
  int n = 0;
  for (const auto& s : outputs)
    if (s.active && s.interval > 0 && std::find(s.fields.begin(), s.fields.end(), field) != s.fields.end()) ++n;
  return n;
}

bool StreamManager::any_ringing(long long now) const {
  for (const auto& s : outputs)
    if (s.ringing(now)) return true;
  return false;
}

const OutputStream* StreamManager::find(const std::string& name) const {
  for (const auto& s : outputs)
    if (s.spec.name == name) return &s;
  return nullptr;
}

OutputStream* StreamManager::find(const std::string& name) {
  for (auto& s : outputs)
    if (s.spec.name == name) return &s;
  return nullptr;
}

StreamManager make_stream_manager(const std::vector<StreamSpec>& specs, const ActivePackages& packages,
                                  const std::vector<std::string>& scalar_names, bool do_restart) {
  StreamManager sm;
  std::set<std::string> names;
  for (const auto& s : specs) {
    if (!names.insert(s.name).second) throw Error("streams: stream '" + s.name + "' defined twice");
    const bool active = any_active(s.packages, packages);
    if (s.input && active) {
      // atm_core_init reads 'input' (a cold start) and atm_core_run reads every input stream whose input alarm
      // rings (mpas_atm_core.F:779-812): only initial_only / none intervals of 'input' and 'restart' are handled
      const bool initial = s.input_interval == "initial_only" || s.input_interval == "none";
      if (!initial || (s.name != "input" && s.name != "restart" && s.input_interval != "none"))
        throw Error("streams: input stream '" + s.name + "' with input_interval '" + s.input_interval +
                    "' has no verified port");
    }
    if (!s.output) continue;
    OutputStream o;
    o.spec = s;
    o.active = active;
    if (s.output_interval != "none") {
      if (s.output_interval == "initial_only")
        throw Error("streams: output stream '" + s.name + "' at 'initial_only' has no verified port");
      o.interval = interval_seconds(s.output_interval);
      if (o.interval <= 0) throw Error("streams: stream '" + s.name + "' has a zero output interval");
    }
    if (!o.active || o.interval == 0) {
      sm.outputs.push_back(std::move(o));
      continue;
    }
    if (s.filename_interval != "" && s.filename_interval != "output_interval")
      throw Error("streams: output stream '" + s.name + "' with filename_interval '" + s.filename_interval +
                  "' has no verified port");
    if (s.immutable && s.name != "restart")
      throw Error("streams: immutable output stream '" + s.name + "' has no verified port");
    if (s.filename_template.empty()) throw Error("streams: output stream '" + s.name + "' has no filename_template");
    std::set<std::string> seen;
    for (const auto& f : s.fields) {
      const RegistryVar* rv = find_registry_var(f);
      std::vector<std::string> expanded;
      if (!rv) {
        // a var_array: its constituents (Registry.xml:1701-1748: the scalars the packages activate)
        if (f != "scalars") throw Error("streams: field '" + f + "' of stream '" + s.name + "' is not in the Registry");
        expanded = scalar_names;
      } else {
        if (!packages.active(*rv)) continue;  // an inactive field is skipped
        expanded.push_back(f);
      }
      for (auto& e : expanded)
        if (seen.insert(e).second) o.fields.push_back(e);
    }
    sm.outputs.push_back(std::move(o));
  }
  // mpas_atm_core.F:667: no restart file at the initial time
  if (OutputStream* r = sm.find("restart")) r->reset(0);
  // mpas_atm_core.F:669-673: a restart run writes no history or diagnostics at its start either (the files of the
  // run it continues). The alarms start at the clock's start, the restart time (reference_time 'initial_time',
  // xml_stream_parser.c:1291-1296), so they first ring one interval after it.
  if (do_restart)
    for (const char* n : {"output", "diagnostics"})
      if (OutputStream* o = sm.find(n)) o->reset(0);
  return sm;
}

}  // namespace cppmpas
