/*
 * ElementPacket.h - Pass around fundamental geometry
 *
 * (c)2018-20 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
 *            Blake B Hillier <blakehillier@mac.com>
 */

#pragma once

#include "Omega2D.h"

#include <algorithm>
#include <vector>

// Helper class for passing arbitrary elements around
template<class S>
class ElementPacket {
public:
  ElementPacket<S>(std::vector<S> _x = std::vector<S>(),
                   std::vector<Int> _idx = std::vector<Int>(),
                   std::vector<S> _val = std::vector<S>(),
                   size_t _nelem = 0,
                   uint8_t _ndim = -1)
    : x(_x), idx(_idx), val(_val), nelem(_nelem), ndim(_ndim)
    {}
  ~ElementPacket<S>() = default;

  ElementPacket<S>(ElementPacket<S> const&) = default; //allow copy
  ElementPacket<S>(ElementPacket<S>&&) = default; //allow move
  ElementPacket<S>& operator=(ElementPacket<S> const&) = default; //allow copy
  ElementPacket<S>& operator=(ElementPacket<S>&&) = default; //allow move

  // Ensure the element packet is correct
  // This will probaly need some parameters and some if statements
  // Check the sim.add_* functions
  bool verify(int attribute, int check) { return (attribute % check == 0); }
  void add(ElementPacket<S> packet) {
    x.insert(x.end(), packet.x.begin(), packet.x.end());
    idx.erase(idx.end()-2, idx.end());
    // Add the last current vertex number to the new set of indices (except the 0 at the end)
    std::transform(packet.idx.begin(), packet.idx.end()-1, packet.idx.begin(),
                   std::bind2nd(std::plus<Int>(), idx.back()));
    idx.insert(idx.end(), packet.idx.begin(), packet.idx.end());
    val.insert(val.end(), packet.val.begin(), packet.val.end());
  }

  std::vector<S> x;
  std::vector<Int> idx;
  std::vector<S> val;
  size_t nelem;
  uint8_t ndim;	// 0=points, 1=surfaces, 2=volumes for 2D
		// 0=points, 1=lines, 2=surfaces, 3=volumes for 3D
};

