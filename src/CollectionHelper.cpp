/*
 * CollectionHelper.cpp - Not sure why I need this, but it won't compile without it, so here you go.
 *
 * (c)2018 Applied Scientific Research, Inc.
 *         Written by Mark J Stock <markjstock@gmail.com>
 */

#include "Collection.h"

#include <iostream>
#include <variant>

// helper function
std::string to_string(Collection& _c) {
  return std::visit([=](auto& elem) { return elem.to_string(); }, _c);
}

