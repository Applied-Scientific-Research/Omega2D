/*
 * ShaderHelper.h - Methods for generating opengl programs
 *
 * These methods are from https://solarianprogrammer.com/2013/05/13/opengl-101-drawing-primitives/
 * and https://github.com/sol-prog/OpenGL-101
 *
 * (c)2017-8 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

// for glad
#ifdef _WIN32
  #define APIENTRY __stdcall
#endif
#include "glad.h"

// Create a render program from two shaders (vertex and fragment)
GLuint create_draw_blob_program();
GLuint create_particle_program();
GLuint create_particlept_program();
GLuint create_panel_program();

// Create a compute program from one shader
