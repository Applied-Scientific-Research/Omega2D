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
#include <functional>
#include <vector>
#include <iostream>

// Helper class for passing arbitrary elements around
template<class S>
class ElementPacket {
public:
  ElementPacket<S>(std::vector<S> _x = std::vector<S>(),
                   std::vector<Int> _idx = std::vector<Int>(),
                   std::vector<S> _val = std::vector<S>(),
                   size_t _nelem = 55,
                   uint8_t _ndim = 55)
    : x(_x), idx(_idx), val(_val), nelem(_nelem), ndim(_ndim)
    {}
  ~ElementPacket<S>() = default;

  ElementPacket<S>(ElementPacket<S> const&) = default; //allow copy
  ElementPacket<S>(ElementPacket<S>&&) = default; //allow move
  ElementPacket<S>& operator=(ElementPacket<S> const&) = default; //allow copy
  ElementPacket<S>& operator=(ElementPacket<S>&&) = default; //allow move

  // Ensure the element packet is correct
  // These functions have been made for boundary segments. They should be tested and abstracted
  // To deal with points, etc. as well
  bool verify(int attribute, int check) {
    bool good = (attribute % check == 0);
    
    if (x.size() % Dimensions != 0) {
      std::cout << "Size of x is not divisible by Dimensions!\nx.size(): ";
      std::cout  << x.size() << " Dimensions: " << Dimensions << std::endl;
      good = false;
    }
    
   /* if (val.size() % x.size()/Dimensions != 0) {
      std::cout << "val does not have correct number of panels!\nval.size(): ";
      std::cout  << val.size() << " number of panels: " << x.size()/Dimensions << std::endl;
      good = false;
    }*/
    return good;
  }
 
  // This simply appends an Element Packet to the end of the existing packet,
  // removing any duplicate points. The end result currently is a longer segment.
  // Boundaries will need to add another set of idx values and a 0 to close it.
  // I don't think this will currently work with flow/measure features, but ideally
  // this would also be able to merge those packets.
  void add(ElementPacket<S> packet) {
    // Check if they have overlapping points on the edges
    int samef = 0;
    for (size_t i = x.size()-1; i > x.size()-Dimensions-1; i--) {
      const int j = (i+Dimensions) % x.size();
      if (x[i] == packet.x[j]) { samef += 1; }
    } // also need to do idx and val. Could just create erase function for packets
    if (samef == Dimensions) { packet.x.erase(packet.x.begin(), packet.x.begin()+Dimensions); }
    int sameb = 0;
    for (size_t i = packet.x.size()-1; i > packet.x.size()-Dimensions-1; i--) {
      const int j = (i+Dimensions) % packet.x.size();
      if (packet.x[i] == x[j]) { sameb += 1; }
    }
    if (sameb == Dimensions) { packet.x.erase(packet.x.end()-Dimensions, packet.x.end()); }
    // Adjust indices if we removed 2 duplicate points
    if ((samef == Dimensions) && (sameb == Dimensions)) {
      packet.idx.erase(packet.idx.end()-Dimensions, packet.idx.end());
    }

    // Combine vectors 
    x.insert(x.end(), packet.x.begin(), packet.x.end());
    // Add the last current vertex number to the new set of indices
    std::transform(packet.idx.begin(), packet.idx.end(), packet.idx.begin(),
                   std::bind(std::plus<Int>(), std::placeholders::_1, idx.back()));
    idx.insert(idx.end(), packet.idx.begin(), packet.idx.end());
    val.insert(val.end(), packet.val.begin(), packet.val.end());
    nelem = val.size();
  }

  std::vector<S> x;
  std::vector<Int> idx;
  std::vector<S> val;
  size_t nelem;
  uint8_t ndim;	// 0=points, 1=surfaces, 2=volumes for 2D
		// 0=points, 1=lines, 2=surfaces, 3=volumes for 3D
};

