/*
 * Omega2D.h - Useful definitions for anywhere in the code
 *
 * (c)2018-20 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
 */

#pragma once

// Use this for indexes into panels or bodies
// this means we can have no more than 65536 element segments in the system - seems reasonable for 2D
#include <cstdint>
using Int = uint16_t;
#include <cstdlib>

const size_t Dimensions = 2;

// element type
enum elem_t {
  active   = 1,  // active vorticity
  reactive = 2,  // active once strength is solved
  inert    = 3,  // does not affect flow
  hybrid   = 4   // does not affect flow
};

// movement type
enum move_t {
  lagrangian = 1, // moves with local velocity
  bodybound  = 2, // moves with attached body
  fixed      = 3  // does not move
};

