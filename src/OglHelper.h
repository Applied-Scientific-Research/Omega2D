/*
 * OglHelper.h - Useful defines for templated OpenGL calls
 *
 * (c)2017-8 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
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

#include "glad.h"

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
