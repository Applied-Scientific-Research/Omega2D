/*
 * ElementPacket.h - Pass around fundamental geometry
 *
 * (c)2018-21 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
 *            Blake B Hillier <blakehillier@mac.com>
 */

#pragma once

#include "Omega2D.h"

#include <algorithm>
#include <functional>
#include <vector>
#include <iostream>
#include <cassert>


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

  // empty packet constructor (forces dimension)
  ElementPacket<S>(uint8_t _ndim)
    : ElementPacket<S>(std::vector<S>(), std::vector<Int>(),
                       std::vector<S>(), 0, _ndim)
    {}

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
  void add(ElementPacket<S> _in) {

    // We shouldn't get here if we're adding elements of different dimensions
    assert(ndim == _in.ndim && "Should not be adding packets of different types");

    // Don't add packets with different numbers of values per element
    assert(val.size()/nelem == _in.val.size()/_in.nelem && "Packet val array size mismatch");

    // No special work required if we're only adding points to points
    //if (ndim == 0) {
    //}

    /* Do no special manipulation - it's OK to have unique, overlapping nodes
    // Check if they have overlapping points on the edges
    int samef = 0;
    for (size_t i = x.size()-1; i > x.size()-Dimensions-1; i--) {
      const size_t j = (i+Dimensions) % x.size();
      if (x[i] == _in.x[j]) { samef += 1; }
    } // also need to do idx and val. Could just create erase function for packets
    if (samef == Dimensions) { _in.x.erase(_in.x.begin(), _in.x.begin()+Dimensions); }

    int sameb = 0;
    for (size_t i = _in.x.size()-1; i > _in.x.size()-Dimensions-1; i--) {
      const int j = (i+Dimensions) % _in.x.size();
      if (_in.x[i] == x[j]) { sameb += 1; }
    }
    if (sameb == Dimensions) { _in.x.erase(_in.x.end()-Dimensions, _in.x.end()); }

    // Adjust indices if we removed 2 duplicate points
    if ((samef == Dimensions) && (sameb == Dimensions)) {
      _in.idx.erase(_in.idx.end()-Dimensions, _in.idx.end());
    }
    */

    // Add the last current vertex number to the new set of indices
    const size_t nodecnt = x.size() / Dimensions;
    std::transform(_in.idx.begin(), _in.idx.end(), _in.idx.begin(),
                   std::bind(std::plus<Int>(), std::placeholders::_1, nodecnt));
    idx.insert(idx.end(), _in.idx.begin(), _in.idx.end());

    // Combine the rest of the vectors
    x.insert(x.end(), _in.x.begin(), _in.x.end());
    val.insert(val.end(), _in.val.begin(), _in.val.end());

    // be careful about this - don't use val.size()
    nelem += _in.nelem;
  }

  void print() {
    for(size_t i=0; i<x.size()/Dimensions; i++) {
      std::cout << "idx: " << idx[i] << " x: " << x[Dimensions*i] << " " << x[Dimensions*i+1] << std::endl;
    }
    std::cout << std::endl;
  }

  std::vector<S> x;
  std::vector<Int> idx;
  std::vector<S> val;
  size_t nelem;
  uint8_t ndim;	// 0=points, 1=surfaces, 2=volumes for 2D
		// 0=points, 1=lines, 2=surfaces, 3=volumes for 3D
};

