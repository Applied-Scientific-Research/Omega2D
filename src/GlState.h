/*
 * GlState.h - Store general OpenGL state for drawing a Collection
 *
 * (c)2019 Applied Scientific Research, Inc.
 *         Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "glad.h"

#include <iostream>
#include <vector>
#include <cassert>


// Hold state for OpenGL
class GlState {
public:
  GlState(const int _nvbo, const int _nspo) {

    assert(_nvbo>=0 && "Invalid number of VBOs requested");
    assert(_nspo>=0 && "Invalid number of shader program objects requested");
    std::cout << "  new GlState with " << _nvbo << " buffers and " << _nspo << " shader programs" << std::endl;

    // generate the vao
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    // and the vbo array
    vbo.resize(_nvbo);
    glGenBuffers(_nvbo, vbo.data());

    // finally the shader program array
    spo.resize(_nspo);

    num_uploaded = 0;
  }

  // must specifically destroy buffers
  ~GlState() {
    //std::cout << "In ~GlState glIsVertexArray(vao) == " << (glIsVertexArray(vao)==GL_TRUE) << std::endl;
    if (glIsVertexArray(vao) == GL_FALSE) return;
    //std::cout << " finishing ~GlState" << std::endl;

    glBindVertexArray(vao);
    glDeleteBuffers(vbo.size(), vbo.data());
    for (size_t i=0; i<spo.size(); ++i) {
      glDeleteProgram(spo[i]);
    }
    glDeleteVertexArrays(1,&vao);
    glBindVertexArray(0);
  }

  //
  // member variables
  //

  // the one vertex array object per Collection
  GLuint vao;

  // some number of vertex buffer objects (usually 4)
  std::vector<GLuint> vbo;
  GLuint qvbo;

  // shader program objects (1 or 2)
  std::vector<GLuint> spo;

  // we probably don't need this
  GLsizei num_uploaded;

  // some number of attributes
  GLint projmat_attribute, projmat_attribute_bl, projmat_attribute_pt, quad_attribute_bl, quad_attribute_pt;
  GLint def_color_attribute, pos_color_attribute, neg_color_attribute, str_scale_attribute, unif_rad_attribute;
  GLint rad_scale_attribute, one_color_attribute;
};

