/*
 * Collection.h - definition of variant type for element collections
 *
 * (c)2018-20 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Points.h"
#include "Surfaces.h"
#include "Volumes.h"

#include <variant>

// alias for any type of collection of elements

using Collection = std::variant<Points<float>,
                                Surfaces<float>,
                                Volumes<float>>;

