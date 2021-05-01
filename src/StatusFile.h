/*
 * StatusFile.h - Class to contain and compose a line-per-step status file
 *
 * (c)2019,21 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include <json/json.hpp>

#include <variant>
#include <string>
#include <vector>

// write a simple space-delimited (gnuplot-friendly) or comma-delimited (Paraview-friendly)
enum StatusFormat { dat, csv };

// a float or integer element
using StatusValue = std::variant<float, int>;


// 1-D elements
class StatusFile {
public:
  StatusFile()
  : use_it(false),
    format(dat),
    num_sims(0),
    num_lines(0),
    fn(""),
    names(),
    vals()
  {}

  // functions
  bool is_active();
  void set_filename(const std::string);
  std::string get_filename();
  void reset_sim();
  void append_value(const std::string,const float);
  void append_value(const std::string,const int);
  void append_value(const float);
  void append_value(const int);
  void write_line();

  // in/out
  void from_json(const nlohmann::json);
  void add_to_json(nlohmann::json&) const;

  // member variables
  bool use_it;
  StatusFormat format;		// format (dat or csv)
  int num_sims;				// number of data sets in this file
  int num_lines;			// number of data lines in this set
  std::string fn;			// the status file name
  std::vector<std::string> names;	// names of the variables
  std::vector<StatusValue> vals;	// values to write at each step
};

