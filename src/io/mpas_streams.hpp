// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// MPAS's stream manager for the output streams of a run, as atm_core_run drives it (mpas_atm_core.F:625-919):
//  * the streams.atmosphere file (xml_stream_parser.c: <immutable_stream> / <stream> with their attributes, the
//    fields of a <stream> from <var>, <var_array> and <file> lists, read from the directory of the streams file) and
//    the immutable streams' fields from the Registry (Registry.xml:602-925, io/registry.hpp);
//  * a stream is active when it has no packages attribute or one of its packages is active
//    (mpas_stream_manager.F:MPAS_stream_mgr_ringing_alarms / stream % active_stream); a field of a stream is active
//    when its own packages are (mpas_stream_manager.F:3430-3480, the field's isActive);
//  * the output alarm of a stream with an output_interval other than 'none' rings at the start time and every
//    interval after it (mpas_stream_manager.F:604-640, reference_time = the start time); a reset moves it past now;
//  * the file name: filename_template with $Y $M $D $h $m $s of the valid time (mpas_timekeeping.F
//    mpas_expand_string), one file per output time (filename_interval = output_interval, the default of an output
//    stream, mpas_stream_manager.F:550-560).
// One rank, one block. Anything else the stream manager can do (input streams other than the initial 'input' and
// 'restart' reads, clobber modes, precision, reference_time, io_type, filename_interval for output streams, output
// at 'initial_only', 'none'-interval streams with active packages that something would write) is an error naming
// the stream and the attribute.
#include <string>
#include <vector>

#include "io/registry_pools.hpp"

namespace cppmpas {

struct StreamSpec {
  std::string name;
  bool immutable = false, input = false, output = false;
  std::string filename_template, filename_interval, input_interval = "none", output_interval = "none", packages;
  std::vector<std::string> fields;  // as listed (var_array names unexpanded); immutable streams: the Registry's
};

// Reads streams.atmosphere (and the <file> lists it names, relative to its directory).
std::vector<StreamSpec> read_streams_file(const std::string& path);

// An MPAS time interval string ("D_hh:mm:ss", "hh:mm:ss", "mm:ss", "ss", optional ".fraction" refused) in seconds.
long long interval_seconds(const std::string& s);

// mpas_expand_string for $Y $M $D $h $m $s; stamp is "YYYY-MM-DD_hh:mm:ss".
std::string expand_filename(const std::string& filename_template, const std::string& stamp);

struct OutputStream {
  StreamSpec spec;
  bool active = false;
  long long interval = 0;  // seconds
  long long ring = 0;      // the next ring time, seconds after the start
  std::vector<std::string> fields;  // active fields, var_arrays expanded to their active constituents
  bool ringing(long long now) const { return active && interval > 0 && now >= ring; }
  void reset(long long now) {  // mpas_reset_clock_alarm: the next ring time after now
    if (interval > 0)
      while (ring <= now) ring += interval;
  }
};

struct StreamManager {
  std::vector<OutputStream> outputs;
  // UT:79-121 MPAS_field_will_be_written at time now
  bool field_will_be_written(const std::string& field, long long now) const;
  // UT:130-213 mpas_stream_inclusion_count(field, direction=MPAS_STREAM_OUTPUT)
  int output_inclusion_count(const std::string& field) const;
  bool any_ringing(long long now) const;
  const OutputStream* find(const std::string& name) const;
  OutputStream* find(const std::string& name);
};

// The output streams of a run: every stream of the file is checked (an input stream other than 'input' and 'restart'
// that is active with an input interval, or an unsupported attribute, is an error); scalar_names: the var_array
// scalars' constituents of the run. The restart stream's output alarm is reset at the start (mpas_atm_core.F:667).
StreamManager make_stream_manager(const std::vector<StreamSpec>& specs, const ActivePackages& packages,
                                  const std::vector<std::string>& scalar_names, bool do_restart);

}  // namespace cppmpas
