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

  // get the boundary corresponding to the wall
  const ReadMsh::boundary wall = mesh.get_bdry("wall");
  if (wall.N_edges == 0) {
    std::cout << "  no boundary called 'wall' in this msh file, skipping." << std::endl;
    return ElementPacket<S>();
  }

  // read *all* the nodes in
  const std::vector<ReadMsh::node>& nodes = mesh.get_nodes();
  for (auto& thisnode : nodes) {
    const ReadMsh::Cmpnts2& thispos = thisnode.coor;
    x.push_back(thispos.x);
    x.push_back(thispos.y);
  }

  // get a reference to the complete edge list
  const std::vector<ReadMsh::edge>& edges = mesh.get_edges();

  // find out how large each array will be
  const size_t np = wall.N_edges;
  std::cout << "wall has " << np << " edges" << std::endl;
  for (uint32_t thisedge : wall.edges) {
    std::cout << "  edge " << thisedge << " has " << edges[thisedge].N_nodes << " nodes" << std::endl;
    // 1st and 2nd nodes are the end nodes, regardless of how many nodes there are on this edge
    assert(edges[thisedge].N_nodes > 1 && "Edge does not have enough nodes!");
    // HACK - annular gmsh meshes have wall defined CCW (right wall is to fluid), not CW (left wall is)
    idx.push_back(edges[thisedge].nodes[1]);
    idx.push_back(edges[thisedge].nodes[0]);
  }

  // compress the nodes vector to remove unused, adjust idx pointers

  std::vector<int32_t> newidx(nodes.size());
  // -1 means that this node is not used
  std::fill(newidx.begin(), newidx.end(), -1);
  // flag all nodes that are used
  for (auto& thisidx : idx) newidx[thisidx] = thisidx;
  // compress the x vector first
  size_t nnodesused = 0;
  for (size_t i=0; i<newidx.size(); ++i) {
    if (newidx[i] == -1) {
      // this node is not used in the wall boundary
    } else {
      // this node *is* used
      // move it backwards
      x[2*nnodesused]   = x[2*i];
      x[2*nnodesused+1] = x[2*i+1];
      // and tell the index where it moved to
      newidx[i] = nnodesused;
      // increment the counter
      nnodesused++;
    }
  }
  x.resize(2*nnodesused);
  // reset indices to indicate their new position in the compressed array
  for (auto& thisidx : idx) thisidx = newidx[thisidx];

  // set boundary condition value to 0.0 (velocity BC)
  vals.resize(np);
  std::fill(vals.begin(), vals.end(), 0.0);

  // return the element packet
  ElementPacket<S> packet({x, idx, vals, (size_t)(np), (uint8_t)1});
  if (packet.verify(packet.x.size(), Dimensions)) {
    return packet;
  } else {
    return ElementPacket<S>();
  }
}
