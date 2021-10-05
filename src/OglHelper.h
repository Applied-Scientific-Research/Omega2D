/*
 * OglHelper.h - Useful defines for templated OpenGL calls
 *
 * (c)2017-8 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
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

// for glad
#ifndef APIENTRY
  #ifdef _WIN32
    #define APIENTRY __stdcall
  #else
    #define APIENTRY
  #endif
  #define GL_APIENTRY_DEFINED
#endif // APIENTRY

#include <glad/glad.h>

#include <cstdint>

template <class T>
struct GetTypeHelper {
};

template <> struct GetTypeHelper<uint8_t> { enum { type = GL_UNSIGNED_BYTE }; };
template <> struct GetTypeHelper<uint16_t> { enum { type = GL_UNSIGNED_SHORT }; };
template <> struct GetTypeHelper<uint32_t> { enum { type = GL_UNSIGNED_INT }; };
template <> struct GetTypeHelper<float> { enum { type = GL_FLOAT }; };
template <> struct GetTypeHelper<double> { enum { type = GL_DOUBLE }; };

template <class T>
constexpr auto get_gl_type = GetTypeHelper<T>::type;
