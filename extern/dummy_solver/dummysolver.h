//
// dummysolver.h
//
// (c)2020 Applied Scientific Research, Inc.
//         Mark J. Stock <markjstock@gmail.com>
#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdint>

namespace DummySolver {

class Solver {

public:
  Solver() {}  //default constructor
  ~Solver() {}  //destructor define here later

  const uint32_t get_nnodes() const { return N_nodes; }
  const uint32_t get_nelems() const { return N_elements; }

  // only allow const getters
  const std::vector<double>& get_nodes() { return nodes; }
  const std::vector<uint32_t>& get_elems() { return elems; }

  // read the given mesh file and populate the data structures
  int32_t read_msh_file(const char* const filename);

  // find a searchname in the MSH file and return the file stream after that line
  //bool locate_in_file(std::ifstream& filestream, const std::string& searchname);

private:
  uint32_t N_nodes = 0;     //# of nodes read from msh file
  uint32_t N_elements = 0;  //# of elements read from msh file
  std::vector<double> nodes;  // coordinates of the nodes
  std::vector<uint32_t> elems;  // coordinates of the nodes

}; // end class Solver

}  // end namespace
