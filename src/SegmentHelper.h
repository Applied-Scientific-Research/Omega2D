/*
 * SegmentHelper.h - Useful routines for segments (group of edges or Surface elems)
 *
 * (c)2021 Applied Scientific Research, Inc.
 *         Written by Mark J Stock <markjstock@gmail.com>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
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

