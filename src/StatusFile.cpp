/*
 * StatusFile.cpp - Class to contain and compose a line-per-step status file
 *
 * (c)2019 Applied Scientific Research, Inc.
 *         Written by Mark J Stock <markjstock@gmail.com>
 */

#include "StatusFile.h"

#include <iostream>
#include <fstream>
//#include <cassert>


// accept a filename, means that we should use it
void
StatusFile::set_filename(const std::string _fn) {
  use_it = true;
  fn = _fn;
}

std::string
StatusFile::get_filename() {
  return fn;
}

// begin writing a new data set to the file
void
StatusFile::reset_sim() {
  if (use_it) {
    num_lines = 0;
  }
}

// append an int or float to the list to write
void
StatusFile::append_value(const float _val) {
  vals.push_back(_val);
}
void
StatusFile::append_value(const int _val) {
  vals.push_back(_val);
}

// write a line
void
StatusFile::write_line() {
  if (use_it) {
    std::ofstream outfile;
    outfile.open(fn, std::ios::app);

    // is this a new data set?
    if (num_lines == 0) {
      if (num_sims > 0) outfile << std::endl;
      outfile << "# new run" << std::endl;
      num_sims++;
    }

    // write data line
    for (auto &val : vals) {
      std::visit([&outfile](const auto& v) { outfile << v; }, val);
      if (val == vals.back()) {
        outfile << std::endl;
      } else {
        outfile << " ";
      }
    }
    num_lines++;
  }

  // empty the vector
  vals.clear();
}

