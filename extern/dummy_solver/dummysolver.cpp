//
// dummysolver.cpp
//
// (c)2020 Applied Scientific Research, Inc.
//         Mark J. Stock <markjstock@gmail.com>

#include "dummysolver.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cassert>
#include <cstdint>

namespace DummySolver {


int32_t Solver::read_msh_file (const char* const filename) {
  // reads a msh file output from the Gmsh software. The msh file is in ASCII 4.1 version of the Gmsh output
  std::cout << "     Gmsh file   ***** " << filename << " *****   opened for reading ..." << std::endl << std::endl;
  int32_t retval = 1, tmp, tmp1, tmp2;
  uint32_t index;
  uint32_t element_type, N_boundary, node_index;
  uint32_t entities_N_points, entities_N_curves, entities_N_faces;
  double double_field;
  std::vector<uint32_t> unorganized_node_mapping; // maps the node numbers read from the msh file into the unorganized indices (helps to remove the indices gap)
  std::vector<uint8_t> unorganized_node_type; //0 for corner, 1 for edge and 2 for nodes on the faces
  std::ifstream mshfile(filename);
  if (mshfile.fail()) {
    std::cout << "Input file opening failed.\n";
    exit(1);
  }

  return retval;
}

} // end namespace
