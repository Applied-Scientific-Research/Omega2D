/*
 * ShaderHelper.h - Methods for generating opengl programs
 *
 * These methods are from https://solarianprogrammer.com/2013/05/13/opengl-101-drawing-primitives/
 * and https://github.com/sol-prog/OpenGL-101
 *
 * (c)2017-20 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
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

#include <string>

// Create a render program from two shaders (vertex and fragment)
GLuint create_draw_blob_program();
GLuint create_draw_point_program();
GLuint create_draw_surface_line_prog();
GLuint create_vertfrag_prog(const std::string, const std::string);
//GLuint create_vertgeomfrag_prog(const std::string, const std::string, const std::string);

// Create a compute program from one shader
