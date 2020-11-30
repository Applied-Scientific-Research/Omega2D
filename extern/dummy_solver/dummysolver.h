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
  void set_re_d(const double);
  void set_elemorder_d(const int32_t);
  void set_timeorder_d(const int32_t);
  void set_numsteps_d(const int32_t);
  void init_d(std::vector<double>, std::vector<uint32_t>, std::vector<uint32_t>, std::vector<uint32_t>);
  std::vector<double> getsolnpts_d();
  std::vector<double> getopenpts_d();
  void setopenvels_d(std::vector<double>);
  void setsolnvort_d(std::vector<double>);
  void solveto_d(const double, const int32_t, const int32_t, const double);
  std::vector<double> getallvorts_d();

private:

  // the physical geometry
  uint32_t N_nodes = 0;     //# of nodes read from msh file
  uint32_t N_elements = 0;  //# of elements read from msh file
  std::vector<double> nodes;  // coordinates of the nodes
  std::vector<uint32_t> elems;  // pointers to nodes for all elements (4 or 9 per elem)
  std::vector<uint32_t> wbdry;  // pointers to nodes of wall boundary elements (2 or 3 per elem)
  std::vector<uint32_t> obdry;  // pointers to nodes of open boundary elements (2 or 3 per elem)

  // the solution points
  uint32_t N_snodes = 0;     //# of nodes in the solution space (order+1)^2 per element
  uint32_t N_selements = 0;  //# of 2d elements in the solution space
  std::vector<double> snodes;  // coordinates of the solution nodes
  std::vector<uint32_t> selems;  // pointers to nodes for all solution elements
  std::vector<uint32_t> sopts;  // pointers to nodes for all solution elements on open boundaries

  int32_t num_substeps;	// number of internel substeps per external time step
  int32_t elem_order;	// internal element order
  int32_t time_order;	// internal time integration order
  double curr_time;	// current simulation time (as far as we know it)

  double reynolds;	// reynolds number

}; // end class Solver

}  // end namespace
