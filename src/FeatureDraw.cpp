/*
 * FeatureDraw.cpp - Draw Features before they are initialized
 *
 * (c)2020 Applied Scientific Research, Inc.
 *         Mark J Stock <markjstock@gmail.com>
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

#include "FeatureDraw.h"
#include "ShaderHelper.h"
#include "OglHelper.h"

#include <iostream>
#include <vector>
#include <cassert>

const std::string vert_shader_source =
#include "shaders/elempack.vert"
;
const std::string frag_shader_source =
#include "shaders/elempack.frag"
;

// Control storage and drawing of features before Simulation takes over

// empty out element packet
void FeatureDraw::clear_elements() {
  m_geom.x.clear();
  m_geom.idx.clear();
  m_geom.val.clear();
  m_idx.clear();
}

float enable(float a) { return a*4; }
float disable(float a) { return a/4; }

// add elements to the local list
void FeatureDraw::add_elements(const ElementPacket<float> _in, const bool _enabled) {
  //std::cout << "  FeatureDraw appending " << (_in.x.size()/2) << " nodes and " << (_in.idx.size()/2) << " elements" << std::endl;

  // remember the number of points in the position array before appending
  const Int np_old = m_geom.x.size() / 2;
  m_geom.x.insert(std::end(m_geom.x), std::begin(_in.x), std::end(_in.x));

  // append indices now
  const Int ni_old = m_geom.idx.size();
  m_geom.idx.insert(std::end(m_geom.idx), std::begin(_in.idx), std::end(_in.idx));

  // but shift the new indices
  for (auto it=m_geom.idx.begin()+ni_old; it!=m_geom.idx.end(); ++it) {
    (*it) += np_old;
  }

  // and set the str/val array for per-element color/strength
  const size_t nv_old = m_geom.val.size();
  m_geom.val.insert(std::end(m_geom.val), std::begin(_in.val), std::end(_in.val));
  if (!_enabled) {
    std::for_each(m_geom.val.begin()+nv_old, m_geom.val.end(), &disable);
  }

  // append to m_idx the start and end indices for this ElementPacket
  m_idx.push_back(std::make_pair((int)nv_old, (int)m_geom.val.size()));
  m_vals_changed = true;
  //std::cout << "  FeatureDraw now has " << (m_geom.x.size()/2) << " nodes and " << (m_geom.idx.size()/2) << " elements" << std::endl;
  //std::cout << "  FeatureDraw pushed pair " << m_idx.back().first << " and " << m_idx.back().second << std::endl;
}


// toggle a set of elements to a new "enabled" setting
/*void FeatureDraw::reset_enabled(const size_t _ifeat, const bool _enabled) {
  assert(_ifeat < m_idx.size() && "Enable feature index larger than number of features");
  assert(m_idx[_ifeat].second <= (int)m_geom.val.size() && "Draw array sizes larger than expected");

  // do the work here
  // The state of _enabled is the new state
  if (_enabled) {
    std::for_each(m_geom.val.begin()+m_idx[_ifeat].first, m_geom.val.begin()+m_idx[_ifeat].second, &enable);
  } else {
    std::cout << "disabling..." << m_idx[_ifeat].first << " " << m_idx[_ifeat].second << std::endl;
    std::for_each(m_geom.val.begin()+m_idx[_ifeat].first, m_geom.val.begin()+m_idx[_ifeat].second, &disable);
  }

  m_vals_changed = true;
}*/


// helper function to clean up initGL
void FeatureDraw::prepare_opengl_buffer(const GLuint _prog, const GLuint _idx,
                                        const GLchar* _name, const GLuint _nper) {
  glBindBuffer(GL_ARRAY_BUFFER, m_gl->vbo[_idx]);
  const GLint attribute_pos = glGetAttribLocation(_prog, _name);
  // Specify how the data for position can be accessed
  glVertexAttribPointer(attribute_pos, _nper, get_gl_type<float>, GL_FALSE, 0, 0);
  // Enable the attribute
  glEnableVertexAttribArray(attribute_pos);
}

void FeatureDraw::initGL(std::vector<float>& _projmat,
                         float*              _poscolor,
                         float*              _negcolor,
                         float*              _defcolor,
                         bool                _oneColor) {

  std::cout << "inside FeatureDraw::initGL" << std::endl;

  // generate the opengl state object and bind the vao
  m_gl = std::make_unique<GlState>(3,1);

  // Allocate space, but don't upload the data from CPU to GPU yet
  glBindBuffer(GL_ARRAY_BUFFER, m_gl->vbo[0]);
  glBufferData(GL_ARRAY_BUFFER, 0, m_geom.x.data(), GL_STATIC_DRAW);

  glBindBuffer(GL_ARRAY_BUFFER, m_gl->vbo[1]);
  glBufferData(GL_ARRAY_BUFFER, 0, m_geom.val.data(), GL_STATIC_DRAW);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_gl->vbo[2]);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, 0, m_geom.idx.data(), GL_STATIC_DRAW);

  // Load and create the line-drawing shader program
  //m_gl->spo[0] = create_draw_surface_line_prog();
  m_gl->spo[0] = create_vertfrag_prog(vert_shader_source, frag_shader_source);

  // Now do the two arrays
  prepare_opengl_buffer(m_gl->spo[0], 0, "pos", 2);
  prepare_opengl_buffer(m_gl->spo[0], 1, "str", 1);

  // Get the location of the attributes that enters in the vertex shader
  m_gl->projmat_attribute = glGetUniformLocation(m_gl->spo[0], "Projection");

  // upload the projection matrix
  glUniformMatrix4fv(m_gl->projmat_attribute, 1, GL_FALSE, _projmat.data());

  // locate where the colors and color scales go
  m_gl->pos_color_attribute = glGetUniformLocation(m_gl->spo[0], "pos_color");
  m_gl->neg_color_attribute = glGetUniformLocation(m_gl->spo[0], "neg_color");
  m_gl->def_color_attribute = glGetUniformLocation(m_gl->spo[0], "def_color");
  m_gl->use_def_attribute = glGetUniformLocation(m_gl->spo[0], "use_def");
  //m_gl->use_back_attribute = glGetUniformLocation(m_gl->spo[0], "use_back");

  // send the current values
  glUniform4fv(m_gl->pos_color_attribute, 1, (const GLfloat *)_poscolor);
  glUniform4fv(m_gl->neg_color_attribute, 1, (const GLfloat *)_negcolor);
  glUniform4fv(m_gl->def_color_attribute, 1, (const GLfloat *)_defcolor);
  glUniform1i(m_gl->use_def_attribute, (const GLint)_oneColor);
  //std::cout << "init pos color as " << _poscolor[0] << " " << _poscolor[1] << " " << _poscolor[2] << " " << _poscolor[3] << std::endl;

  // and indicate the fragment color output
  glBindFragDataLocation(m_gl->spo[0], 0, "frag_color");

  glBindVertexArray(0);
}

void FeatureDraw::updateGL() {

  //std::cout << "inside FeatureDraw::updateGL" << std::endl;

  // has this been init'd yet?
  if (not m_gl) return;
  if (glIsVertexArray(m_gl->vao) == GL_FALSE) return;

  // has the number of elements changed?
  //if (m_vals_changed) {//if (m_gl->num_uploaded != (GLsizei)m_geom.idx.size()) {
    glBindVertexArray(m_gl->vao);

    // Indicate and upload the data from CPU to GPU

    // the positions
    glBindBuffer(GL_ARRAY_BUFFER, m_gl->vbo[0]);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float)*m_geom.x.size(), m_geom.x.data(), GL_DYNAMIC_DRAW);

    // strengths
    glBindBuffer(GL_ARRAY_BUFFER, m_gl->vbo[1]);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float)*m_geom.val.size(), m_geom.val.data(), GL_DYNAMIC_DRAW);

    // and element indices
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_gl->vbo[2]);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(Int)*m_geom.idx.size(), m_geom.idx.data(), GL_DYNAMIC_DRAW);
    glBindVertexArray(0);

    // must tell draw call how many elements are there - or, really, how many indices
    m_gl->num_uploaded = m_geom.idx.size();

  /*} else if (m_vals_changed) {
    // just upload the new vals, nothing else
    glBindVertexArray(m_gl->vao);

    // strengths
    glBindBuffer(GL_ARRAY_BUFFER, m_gl->vbo[1]);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float)*m_geom.val.size(), m_geom.val.data(), GL_DYNAMIC_DRAW);

    glBindVertexArray(0);
  }*/
}

void FeatureDraw::drawGL(std::vector<float>& _projmat,
                         RenderParams&       _rparams,
                         bool                _oneColor) {

  //std::cout << "inside FeatureDraw::drawGL" << std::endl;

  // has this been init'd yet?
  if (not m_gl) {
    initGL(_projmat, _rparams.pos_circ_color,
                     _rparams.neg_circ_color,
                     _rparams.default_color,
                     _oneColor);
    updateGL();
  }

  if (/*(m_gl->num_uploaded != (GLsizei)m_geom.idx.size()) ||*/ m_vals_changed) {
    updateGL();
  }

  // draw if there are any elems to draw
  if (m_gl->num_uploaded > 0) {
    glBindVertexArray(m_gl->vao);

    // get blending ready
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ONE);

    // draw as lines
    glUseProgram(m_gl->spo[0]);

    // upload the current projection matrix
    glUniformMatrix4fv(m_gl->projmat_attribute, 1, GL_FALSE, _projmat.data());

    // upload the current color values
    glUniform4fv(m_gl->pos_color_attribute, 1, (const GLfloat *)_rparams.pos_circ_color);
    glUniform4fv(m_gl->neg_color_attribute, 1, (const GLfloat *)_rparams.neg_circ_color);
    glUniform4fv(m_gl->def_color_attribute, 1, (const GLfloat *)_rparams.default_color);
    glUniform1i(m_gl->use_def_attribute, (const GLint)_oneColor);

    // the one draw call here
    glDrawElements(GL_LINES, m_gl->num_uploaded, get_gl_type<Int>, 0);

    // return state
    glEnable(GL_DEPTH_TEST);
    glDisable(GL_BLEND);

    glBindVertexArray(0);
  }
}

