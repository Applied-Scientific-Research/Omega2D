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

  // functions that an external driver program should call
  void set_re_d_(const double _re) { re = _re; }
  void init_d_(std::vector<double>, std::vector<uint32_t>, std::vector<uint32_t>, std::vector<uint32_t>);
  std::vector<double> getsolnpts_d_();
  std::vector<double> getopenpts_d_();
  void setopenvels_d_(std::vector<double>);
  void solveto_d_(const double);
  std::vector<double> getallvorts_d_();

private:
  uint32_t N_nodes = 0;     //# of nodes read from msh file
  uint32_t N_elements = 0;  //# of elements read from msh file
  std::vector<double> nodes;  // coordinates of the nodes
  std::vector<uint32_t> elems;  // coordinates of the nodes

  double re;		// reynolds number

}; // end class Solver

}  // end namespace
