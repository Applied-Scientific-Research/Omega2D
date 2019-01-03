/*
 * Collection.h - definition of variant type for element collections
 *
 * (c)2018 Applied Scientific Research, Inc.
 *         Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Points.h"
//#include "Panels.h"

#include <variant>

// alias for any type of collection of elements
using Collection = std::variant<Points<float>>;

//using Collection = std::variant<Points<float>,
//				Panels<float>>;

