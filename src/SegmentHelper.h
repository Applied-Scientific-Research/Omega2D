/*
 * SegmentHelper.h - Useful routines for segments (group of edges or Surface elems)
 *
 * (c)2021 Applied Scientific Research, Inc.
 *         Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include <vector>

//
// given a set of surface elements (edges), put them in order end-to-end
//   and find the flow through each element based on its parametric position
//
std::vector<float> find_parabolic_vels(const std::vector<float>&,
                                       const std::vector<Int>&,
                                       const size_t);

