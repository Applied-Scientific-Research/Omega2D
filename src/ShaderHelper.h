/*
 * ShaderHelper.h - Methods for generating opengl programs
 *
 * These methods are from https://solarianprogrammer.com/2013/05/13/opengl-101-drawing-primitives/
 * and https://github.com/sol-prog/OpenGL-101
 *
 * (c)2017-20 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
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

#include <string>

// Create a render program from two shaders (vertex and fragment)
GLuint create_draw_blob_program();
GLuint create_draw_point_program();
GLuint create_draw_surface_line_prog();
GLuint create_vertfrag_prog(const std::string, const std::string);
//GLuint create_vertgeomfrag_prog(const std::string, const std::string, const std::string);

// Create a compute program from one shader
