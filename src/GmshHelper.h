/*
 * GmshHelper.h - Non-class helper functions to convert gmsh data
 *
 * (c)2020 Applied Scientific Research, Inc.
 *         Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Omega2D.h"
#include "ElementPacket.h"
#include "read_MSH_Mesh.h"

#include <string>
#include <cstdint>

template <class S>
ElementPacket<S> get_wall_bdry(const std::string _infile) {

  // read gmsh file
  ReadMsh::Mesh mesh;
  int32_t retval = mesh.read_msh_file(_infile.c_str());
  std::cout << "MSH file (" << _infile << ")";
  if (retval == 1) {
    std::cout << " contains " << mesh.get_nnodes() << " nodes";
    std::cout << " and " << mesh.get_nelems() << " elems" << std::endl;
  } else {
    std::cout << " does not exist or did not read properly (code=";
    std::cout << retval << "), skipping." << std::endl;
    return ElementPacket<S>();
  }

  // prepare the data arrays for the element packet
  std::vector<S> x;
  std::vector<Int> idx;
  std::vector<S> vals;

  // get a reference to all node locations
  const std::vector<ReadMsh::node>& node = mesh.get_nodes();
  const std::vector<ReadMsh::edge>& edge = mesh.get_edges()

  // iterate through the boundary elements
  const ReadMsh::boundary wall = mesh.get_boundary("wall");
  if (wall.N_edges == 0) {
    std::cout << "  no boundary called 'wall' in this msh file, skipping." << std::endl;
    return ElementPacket<S>();
  }

  // find out how large each array will be
  const size_t np = wall.N_edges;

  //x.resize(Dimensions*nnt*nnr);
  //const size_t nodesperelem = std::pow(m_order+1, 2);
  //idx.resize(nodesperelem*m_nt*m_nr);

  // pull data from the mesh object

  // return the element packet
  ElementPacket<S> packet({x, idx, vals, (size_t)(np), (uint8_t)1});
  if (packet.verify(packet.x.size(), Dimensions)) {
    return packet;
  } else {
    return ElementPacket<S>();
  }
}
