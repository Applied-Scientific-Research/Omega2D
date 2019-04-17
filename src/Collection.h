/*
 * Collection.h - definition of variant type for element collections
 *
 * (c)2018-9 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Points.h"
#include "Surfaces.h"

#include <variant>

// alias for any type of collection of elements
// eventually will have Volumes<float> here and in 3D
//                  and Lines<float> in 3D only

using Collection = std::variant<Points<float>,
                                Surfaces<float>>;

