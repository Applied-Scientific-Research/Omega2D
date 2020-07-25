/*
 * FeatureDraw.h - Draw Features before they are initialized
 *
 * (c)2020 Applied Scientific Research, Inc.
 *         Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "ElementPacket.h"
#include "RenderParams.h"
#include "GlState.h"


// Control storage and drawing of features before Simulation takes over
class FeatureDraw {

public:
  FeatureDraw() : m_vals_changed(false) { }

  // deleting GlState will destroy buffers
  ~FeatureDraw() { }

  // visible member functions
  void clear_elements();
  void add_elements(const ElementPacket<float>, const bool);
  void reset_enabled(const size_t, const bool);

  void updateGL();
  void drawGL(std::vector<float>&, RenderParams&, bool);

private:
  //
  // member variables
  //
  void prepare_opengl_buffer(const GLuint, const GLuint, const GLchar*, const GLuint);
  void initGL(std::vector<float>&, float*, float*, float*, bool);

  // Collected geometry
  ElementPacket<float> m_geom;
  std::vector<std::pair<int,int>> m_idx;
  bool m_vals_changed;

  // VAO, VBOs, etc.
  std::unique_ptr<GlState> m_gl;
};

