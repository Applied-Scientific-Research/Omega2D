/*
 * Points.h - Specialized class for 2D points
 *
 * (c)2018-20 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Omega2D.h"
#include "VectorHelper.h"
#include "ElementBase.h"
#include "VtkXmlWriter.h"

#ifdef USE_GL
#include "GlState.h"
#include "RenderParams.h"
#include "OglHelper.h"
#include "ShaderHelper.h"
#include <glad/glad.h>
#endif

#include <iostream>
#include <iomanip> // for setfill and setw
#include <vector>
#include <array>
#include <memory>
#include <optional>
#include <random>
#include <cassert>
#include <cmath>
#include <algorithm>


// 0-D elements
template <class S>
class Points: public ElementBase<S> {
public:
  // only constructor now - using ElementPacket
  Points(const ElementPacket<S>& _in,
         const elem_t _e,
         const move_t _m,
         std::shared_ptr<Body> _bp,
         const float _vd)
    : ElementBase<S>(0, _e, _m, _bp),
      max_strength(-1.0) {

    // ensure that this packet really is Points
    assert(_in.idx.size() == 0 && "Input ElementPacket is not Points");
    assert(_in.ndim == 0 && "Input ElementPacket is not Points");

    // and that it has the right number of values per particle
    if (_e == inert) assert(_in.val.size() == 0 && "Input ElementPacket with fldpts has val array");
    else if (_e == reactive) assert(false && "Input ElementPacket with reactive points is unsupported");
    else assert(_in.val.size() == _in.nelem && "Input ElementPacket with vortons has incorrect size val array");

    // tell the world that we're legit
    std::cout << "  new collection with " << (_in.nelem);
    std::cout << ((_e == inert) ? " tracer" : " vortex") << " elems" << std::endl;
    //std::cout << "  contains " << std::endl;
    //for (size_t i=0; i<_in.nelem; ++i) {
    //  std::cout << "    " << _in.x[2*i] << " " << _in.x[2*i+1] << " " << _in.val[2*i] << " " << _in.val[2*i+1] << std::endl;
    //}

    // need to reset the base class n, because this gets run before the base ctor
    this->n = _in.nelem;

    // this initialization specific to Points - is it, though?
    for (size_t d=0; d<Dimensions; ++d) {
      this->x[d].resize(this->n);
      for (size_t i=0; i<this->n; ++i) {
        this->x[d][i] = _in.x[Dimensions*i+d];
      }
    }

    // save untransformed positions if we are given a Body pointer
    if (_bp) {
      this->ux = this->x;
    }

    if (_e == inert) {
      // field points need neither radius nor strength

    } else {
      // active vortons need radius
      r.resize(this->n);
      std::fill(r.begin(), r.end(), _vd);

      // optional strength in base class
      // need to assign it a vector first!
      // This needs to be modeled after 3D code
      Vector<S> new_s;
      new_s.resize(this->n);
      std::copy(_in.val.begin(), _in.val.end(), new_s.begin());
      this->s = std::move(new_s);
    }

    // velocity in base class
    for (size_t d=0; d<Dimensions; ++d) {
      this->u[d].resize(this->n);
    }
  }

  const Vector<S>& get_rad() const { return r; }
  Vector<S>&       get_rad()       { return r; }

  const S get_averaged_max_str() const { return max_strength; }

  // a little logic to see if we should augment the BEM equations for this object (see Surfaces.h)
  const bool is_augmented() const { return false; }

  // find out the next row index in the BEM after this collection
  // once we start supporting BEM unknowns on points, we'll have to change these
  void set_first_row(const Int _i) { return; }
  const Int get_first_row() const { return 0; }
  const Int get_num_rows()  const { return 0; }
  const Int get_next_row()  const { return 0; }

  const float get_max_bc_value() const { return 0.0; }

  // append more elements this collection
  void add_new(const ElementPacket<S>& _in, const float _vd) {
    // ensure that this packet really is Points
    assert(_in.idx.size() == 0 && "Input ElementPacket is not Points");
    assert(_in.ndim == 0 && "Input ElementPacket is not Points");

    // and that it has the right number of values per particle
    if (this->E == inert) assert(_in.val.size() == 0 && "Input ElementPacket with fldpts has val array");
    else if (this->E == reactive) assert("Input ElementPacket with reactive points is unsupported");
    else assert(_in.val.size() == _in.nelem && "Input ElementPacket with vortons has bad sized val array");

    // remember old size and incoming size (note that Points nelems = nnodes)
    const size_t nold = this->n;
    const size_t nnew = _in.nelem;
    std::cout << "  adding " << nnew << " particles to collection..." << std::endl;

    // must explicitly call the method in the base class first - this pulls out positions and strengths
    ElementBase<S>::add_new(_in);

    // then do local stuff
    if (this->E == inert) {
      // no radius needed

    } else {
      r.resize(nold+nnew);
      std::fill(r.begin()+nold, r.end(), _vd);
    }

    // save the new untransformed positions if we have a Body pointer
    if (this->B) {
      for (size_t d=0; d<Dimensions; ++d) {
        (*this->ux)[d].resize(nold+nnew);
        for (size_t i=nold; i<nold+nnew; ++i) {
          (*this->ux)[d][i] = this->x[d][i];
        }
      }
    }
  }

  // up-size all arrays to the new size, filling with sane values
  void resize(const size_t _nnew) {
    const size_t currn = this->n;
    //std::cout << "  inside Points::resize with " << currn << " " << _nnew << std::endl;

    // must explicitly call the method in the base class - this sets n
    ElementBase<S>::resize(_nnew);

    if (_nnew == currn) return;

    // radii here
    if (this->E == inert) {
      // no radii
    } else {
      const size_t thisrn = r.size();
      r.resize(_nnew);
      for (size_t i=thisrn; i<_nnew; ++i) {
        r[i] = 1.0;
      }
    }
  }

  void zero_vels() {
    // must explicitly call the method in the base class to zero the vels
    ElementBase<S>::zero_vels();
  }

  void finalize_vels(const std::array<double,Dimensions>& _fs) {
    // must explicitly call the method in the base class, too
    ElementBase<S>::finalize_vels(_fs);
  }

  void transform(const double _time) {
    // must explicitly call the method in the base class
    ElementBase<S>::transform(_time);
    // no other specialization required here
  }

  //
  // 1st order Euler advection and stretch
  //
  void move(const double _time, const double _dt) {
    // must explicitly call the method in the base class
    ElementBase<S>::move(_time, _dt);

    // and update the max strength measure
    (void) update_max_str();
  }

  //
  // 2nd order RK advection and stretch
  //
  void move(const double _time, const double _dt,
            const double _wt1, Points<S> const & _u1,
            const double _wt2, Points<S> const & _u2) {
    // must explicitly call the method in the base class
    ElementBase<S>::move(_time, _dt, _wt1, _u1, _wt2, _u2);

    // must confirm that incoming time derivates include velocity (?)

    // and update the max strength measure
    (void) update_max_str();
  }

  // find the new peak strength magnitude
  void update_max_str() {
    S thismax = ElementBase<S>::get_max_str();

    // and slowly update the current saved value
    if (max_strength < 0.0) {
      max_strength = thismax;
    } else {
      max_strength = 0.05*thismax + 0.95*max_strength;
    }
  }

  // add and return the total impulse of all elements
  std::array<S,Dimensions> get_total_impulse() {

    // here is the return vector
    std::array<S,Dimensions> imp;
    imp.fill(0.0);

    if (this->s) {
      // accumulate impulse from each particle
      for (size_t i=0; i<this->n; ++i) {
        imp[0] -= (*this->s)[i] * this->x[1][i];
        imp[1] += (*this->s)[i] * this->x[0][i];
      }
    }

    return imp;
  }

  void add_body_motion(const S _factor, const double _time) {
    // no need to call base class now
    //ElementBase<S>::add_body_motion(_factor);
    // apply a factor times the body motion
    //for (size_t i=0; i<this->get_n(); ++i) {
    //  std::array<S,Dimensions> thisvel = B->get_vel(_time);
    //  // apply the velocity
    //  for (size_t d=0; d<Dimensions; ++d) {
    //    this->u[d][i] += _factor * thisvel[d];
    //  }
    //}
  }

  // augment the strengths with a value equal to that which accounts for
  //   the solid-body rotation of the object
  // do nothing here, but when we get Kutta points, we may need to
  void add_unit_rot_strengths() {}
  void add_solved_rot_strengths(const S _factor) {}

#ifdef USE_GL
  //
  // OpenGL functions
  //

  // helper function to clean up initGL
  void prepare_opengl_buffer(GLuint _prog, GLuint _idx, const GLchar* _name) {
    glBindBuffer(GL_ARRAY_BUFFER, mgl->vbo[_idx]);
    const GLint position_attribute = glGetAttribLocation(_prog, _name);
    // Specify how the data for position can be accessed
    glVertexAttribPointer(position_attribute, 1, get_gl_type<S>, GL_FALSE, 0, 0);
    // Enable the attribute
    glEnableVertexAttribArray(position_attribute);
    // and tell it to advance two primitives per point (2 tris per quad)
    glVertexAttribDivisor(position_attribute, 1);
  }

  // this gets done once - load the shaders, set up the vao
  void initGL(std::vector<float>& _projmat,
              float*              _poscolor,
              float*              _negcolor,
              float*              _defcolor) {

    //std::cout << "inside Points.initGL" << std::endl;
    std::cout << "inside Points.initGL with E=" << this->E << " and M=" << this->M << std::endl;

    // generate the opengl state object with space for 4 vbos and 2 shader programs
    mgl = std::make_shared<GlState>(4,2);

    // Allocate space, but don't upload the data from CPU to GPU yet
    for (size_t i=0; i<Dimensions; ++i) {
      glBindBuffer(GL_ARRAY_BUFFER, mgl->vbo[i]);
      glBufferData(GL_ARRAY_BUFFER, 0, this->x[i].data(), GL_STATIC_DRAW);
    }

    // here is where we split on element type: active/reactive vs. inert
    if (this->E == inert) {

      // Load and create the blob-drawing shader program
      mgl->spo[0] = create_draw_point_program();

      // Only send position arrays - no radius or strength
      prepare_opengl_buffer(mgl->spo[0], 0, "px");
      prepare_opengl_buffer(mgl->spo[0], 1, "py");

      // and for the compute shaders!

      // Get the location of the attributes that enters in the vertex shader
      mgl->projmat_attribute_pt = glGetUniformLocation(mgl->spo[0], "Projection");

      // upload the projection matrix
      glUniformMatrix4fv(mgl->projmat_attribute_pt, 1, GL_FALSE, _projmat.data());

      // locate where the colors and color scales go
      mgl->def_color_attribute = glGetUniformLocation(mgl->spo[0], "def_color");

      // send the current values
      glUniform4fv(mgl->def_color_attribute, 1, _defcolor);

      // locate where the point radius goes
      mgl->unif_rad_attribute = glGetUniformLocation(mgl->spo[0], "rad");

      // send the current values
      glUniform1f(mgl->unif_rad_attribute, (const GLfloat)0.01);

      // and indicate the fragment color output
      glBindFragDataLocation(mgl->spo[0], 0, "frag_color");

      // Initialize the quad attributes
      std::vector<float> quadverts = {-1,-1, 1,-1, 1,1, -1,1};
      glGenBuffers(1, &(mgl->qvbo));
      glBindBuffer(GL_ARRAY_BUFFER, mgl->qvbo);
      glBufferData(GL_ARRAY_BUFFER, quadverts.size()*sizeof(float), quadverts.data(), GL_STATIC_DRAW);

      mgl->quad_attribute_pt = glGetAttribLocation(mgl->spo[0], "quad_attr");
      glVertexAttribPointer(mgl->quad_attribute_pt, 2, GL_FLOAT, GL_FALSE, 0, 0);
      glEnableVertexAttribArray(mgl->quad_attribute_pt);

    } else { // this->E is active or reactive

      glBindBuffer(GL_ARRAY_BUFFER, mgl->vbo[2]);
      glBufferData(GL_ARRAY_BUFFER, 0, r.data(), GL_STATIC_DRAW);
      if (this->s) {
        glBindBuffer(GL_ARRAY_BUFFER, mgl->vbo[3]);
        glBufferData(GL_ARRAY_BUFFER, 0, (*this->s).data(), GL_STATIC_DRAW);
      }

      // Load and create the blob-drawing shader program
      mgl->spo[1] = create_draw_blob_program();

      // Now do the four arrays
      prepare_opengl_buffer(mgl->spo[1], 0, "px");
      prepare_opengl_buffer(mgl->spo[1], 1, "py");
      prepare_opengl_buffer(mgl->spo[1], 2, "rad");
      prepare_opengl_buffer(mgl->spo[1], 3, "sx");

      // and for the compute shaders!

      // Get the location of the attributes that enters in the vertex shader
      mgl->projmat_attribute_bl = glGetUniformLocation(mgl->spo[1], "Projection");

      // upload the projection matrix
      glUniformMatrix4fv(mgl->projmat_attribute_bl, 1, GL_FALSE, _projmat.data());

      // locate where the colors and color scales go
      mgl->pos_color_attribute = glGetUniformLocation(mgl->spo[1], "pos_color");
      mgl->neg_color_attribute = glGetUniformLocation(mgl->spo[1], "neg_color");
      mgl->str_scale_attribute = glGetUniformLocation(mgl->spo[1], "str_scale");
      mgl->rad_scale_attribute = glGetUniformLocation(mgl->spo[1], "rad_scale");

      // send the current values
      glUniform4fv(mgl->pos_color_attribute, 1, (const GLfloat *)_poscolor);
      glUniform4fv(mgl->neg_color_attribute, 1, (const GLfloat *)_negcolor);
      glUniform1f (mgl->str_scale_attribute, (const GLfloat)1.0);
      glUniform1f (mgl->rad_scale_attribute, (const GLfloat)1.0);

      // and indicate the fragment color output
      glBindFragDataLocation(mgl->spo[1], 0, "frag_color");

      // Initialize the quad attributes
      std::vector<float> quadverts = {-1,-1, 1,-1, 1,1, -1,1};
      glGenBuffers(1, &(mgl->qvbo));
      glBindBuffer(GL_ARRAY_BUFFER, mgl->qvbo);
      glBufferData(GL_ARRAY_BUFFER, quadverts.size()*sizeof(float), quadverts.data(), GL_STATIC_DRAW);

      mgl->quad_attribute_bl = glGetAttribLocation(mgl->spo[1], "quad_attr");
      glVertexAttribPointer(mgl->quad_attribute_bl, 2, GL_FLOAT, GL_FALSE, 0, 0);
      glEnableVertexAttribArray(mgl->quad_attribute_bl);

    } // end this->E is active or reactive

    glBindVertexArray(0);
  }

  // this gets done every time we change the size of the positions array
  void updateGL() {
    //std::cout << "inside Points.updateGL" << std::endl;

    // has this been init'd yet?
    if (not mgl) return;
    if (glIsVertexArray(mgl->vao) == GL_FALSE) return;

    const size_t vlen = this->x[0].size()*sizeof(S);
    if (vlen > 0) {
      glBindVertexArray(mgl->vao);

      // Indicate and upload the data from CPU to GPU
      for (size_t i=0; i<Dimensions; ++i) {
        // the positions
        glBindBuffer(GL_ARRAY_BUFFER, mgl->vbo[i]);
        glBufferData(GL_ARRAY_BUFFER, vlen, this->x[i].data(), GL_DYNAMIC_DRAW);
      }

      // here is where we split on element type: active/reactive vs. inert
      if (this->E == inert) {

        // no strengths or radii needed or present

      } else { // this->E is active or reactive

        // the strengths
        if (this->s) {
          glBindBuffer(GL_ARRAY_BUFFER, mgl->vbo[3]);
          glBufferData(GL_ARRAY_BUFFER, vlen, (*this->s).data(), GL_DYNAMIC_DRAW);
        }
        // the radii
        glBindBuffer(GL_ARRAY_BUFFER, mgl->vbo[2]);
        glBufferData(GL_ARRAY_BUFFER, vlen, r.data(), GL_DYNAMIC_DRAW);
      }

      glBindVertexArray(0);

      // must tell draw call how many elements are there
      mgl->num_uploaded = this->x[0].size();
    }
  }

  // OpenGL3 stuff to display points, called once per frame
  void drawGL(std::vector<float>& _projmat,
              RenderParams&       _rparams,
              const float         _vdelta) {

    //std::cout << "inside Points.drawGL" << std::endl;

    // has this been init'd yet?
    if (not mgl) {
      initGL(_projmat, _rparams.pos_circ_color,
                       _rparams.neg_circ_color,
                       _rparams.default_color);
      updateGL();
    }

    if (mgl->num_uploaded > 0) {
      glBindVertexArray(mgl->vao);

      // get blending ready
      glDisable(GL_DEPTH_TEST);
      glEnable(GL_BLEND);
      glBlendFunc(GL_ONE, GL_ONE);

      // here is where we split on element type: active/reactive vs. inert
      if (this->E == inert) {

        // draw as small dots
        glUseProgram(mgl->spo[0]);

        glEnableVertexAttribArray(mgl->quad_attribute_pt);

        // upload the current uniforms
        glUniformMatrix4fv(mgl->projmat_attribute_pt, 1, GL_FALSE, _projmat.data());
        glUniform4fv(mgl->def_color_attribute, 1, (const GLfloat *)_rparams.default_color);
        glUniform1f (mgl->unif_rad_attribute, (const GLfloat)(2.5f*_rparams.tracer_size));

        // the one draw call here
        glDrawArraysInstanced(GL_TRIANGLE_FAN, 0, 4, mgl->num_uploaded);

      } else { // this->E is active or reactive

        // draw as colored clouds
        glUseProgram(mgl->spo[1]);

        glEnableVertexAttribArray(mgl->quad_attribute_bl);

        // upload the current projection matrix
        glUniformMatrix4fv(mgl->projmat_attribute_bl, 1, GL_FALSE, _projmat.data());

        // upload the current color values
        const float color_scaling = _rparams.circ_density * std::pow(_vdelta/_rparams.vorton_scale,2) / max_strength;
        glUniform4fv(mgl->pos_color_attribute, 1, (const GLfloat *)_rparams.pos_circ_color);
        glUniform4fv(mgl->neg_color_attribute, 1, (const GLfloat *)_rparams.neg_circ_color);
        glUniform1f (mgl->str_scale_attribute, (const GLfloat)color_scaling);
        glUniform1f (mgl->rad_scale_attribute, (const GLfloat)_rparams.vorton_scale);

        // the one draw call here
        glDrawArraysInstanced(GL_TRIANGLE_FAN, 0, 4, mgl->num_uploaded);
      }

      // return state
      glEnable(GL_DEPTH_TEST);
      glDisable(GL_BLEND);
      glBindVertexArray(0);
    }
  }
#endif

  std::string to_string() const {
    std::string retstr = " " + std::to_string(this->n) + ElementBase<S>::to_string() + " Points";
    return retstr;
  }

  std::string write_vtk(const size_t _index, const size_t _frameno, const double _time) {
    assert(this->n > 0 && "Inside write_vtk with no points");
  
    const bool asbase64 = true;

    bool has_radii = true;
    bool has_strengths = true;
    std::string prefix = "part_";
    if (this->E==inert) {
      has_strengths = false;
      has_radii = false;
      prefix = "fldpt_";
    }
  
    // generate file name
    std::stringstream vtkfn;
    vtkfn << prefix << std::setfill('0') << std::setw(2) << _index << "_" << std::setw(5) << _frameno << ".vtu";
    VtkXmlWriter ptsWriter = VtkXmlWriter(vtkfn.str(), asbase64);
  
    // include simulation time here
    ptsWriter.addElement("FieldData");
    {
      std::map<std::string, std::string> attribs = {{"type",           "Float64"},
                                                    {"Name",           "TimeValue"},
                                                    {"NumberOfTuples", "1"}};
      ptsWriter.addElement("DataArray", attribs);
      Vector<double> time_vec = {_time};
      ptsWriter.writeDataArray(time_vec);
      // DataArray
      ptsWriter.closeElement();
    }
    // FieldData
    ptsWriter.closeElement();
  
    {
      std::map<std::string, std::string> attribs = {{"NumberOfPoints", std::to_string(this->n).c_str()},
                                                    {"NumberOfCells", std::to_string(this->n).c_str()}};
      ptsWriter.addElement("Piece", attribs);
    }
    
    ptsWriter.addElement("Points");
    {
      std::map<std::string, std::string> attribs = {{"NumberOfComponents", "3"},
                                                    {"Name",               "position"},
                                                    {"type",               "Float32"}};
      ptsWriter.addElement("DataArray", attribs);
      Vector<float> pos = ptsWriter.unpackArray(this->x);
      ptsWriter.writeDataArray(pos);
      // DataArray
      ptsWriter.closeElement();
    }
    // Points
    ptsWriter.closeElement();
  
    ptsWriter.addElement("Cells");
    
    // https://discourse.paraview.org/t/cannot-open-vtu-files-with-paraview-5-8/3759
    // apparently the Vtk format documents indicate that connectivities and offsets
    //   must be in Int32, not UIntAnything. Okay...
    {
      std::map<std::string, std::string> attribs = {{"Name", "connectivity"},
                                                    {"type", "Int32"}};
      ptsWriter.addElement("DataArray", attribs);
      Vector<int32_t> v(this->n);
      std::iota(v.begin(), v.end(), 0);
      ptsWriter.writeDataArray(v);
      // DataArray
      ptsWriter.closeElement();
    }
  
    {
      std::map<std::string, std::string> attribs = {{"Name", "offsets"},
                                                    {"type", "Int32"}};
      ptsWriter.addElement("DataArray", attribs);
      Vector<int32_t> v(this->n);
      std::iota(v.begin(), v.end(), 1);
      ptsWriter.writeDataArray(v);
      // DataArray
      ptsWriter.closeElement();
    }
  
    // except these, they can be chars
    {
      std::map<std::string, std::string> attribs = {{"Name", "types"},
                                                    {"type", "UInt8"}};
      ptsWriter.addElement("DataArray", attribs);
      Vector<uint8_t> v(this->n);
      std::fill(v.begin(), v.end(), 1);
      ptsWriter.writeDataArray(v);
      // DataArray
      ptsWriter.closeElement();
    }
    // Cells
    ptsWriter.closeElement();
  
    {
      std::map<std::string, std::string> attribs = {{"Vectors", "velocity"}};

      std::string scalar_list;
      if (has_strengths) scalar_list.append("circulation,");
      if (has_radii) scalar_list.append("radius,");
      if (this->has_vort()) scalar_list.append("vorticity,");
      if (scalar_list.size()>1) {
        scalar_list.pop_back();
        attribs.insert({"Scalars", scalar_list});
      }

      ptsWriter.addElement("PointData", attribs);
    }
  
    if (has_strengths) {
      std::map<std::string, std::string> attribs = {{"Name", "circulation"},
                                                    {"type", "Float32"}};
      ptsWriter.addElement("DataArray", attribs);
      ptsWriter.writeDataArray(*(this->s));
      ptsWriter.closeElement(); // DataArray
    }

    if (has_radii) {
      std::map<std::string, std::string> attribs = {{"Name", "radius"},
                                                    {"type", "Float32"}};
      ptsWriter.addElement("DataArray", attribs);
      ptsWriter.writeDataArray(this->r);
      ptsWriter.closeElement(); // DataArray
    }

    if (this->has_vort()) {
      std::map<std::string, std::string> attribs = {{"Name", "vorticity"},
                                                    {"type", "Float32"}};
      ptsWriter.addElement("DataArray", attribs);
      ptsWriter.writeDataArray(*(this->w));
      ptsWriter.closeElement(); // DataArray
    }

    {
      std::map<std::string, std::string> attribs = {{"NumberOfComponents", "3"},
                                                    {"Name",               "velocity"},
                                                    {"type",               "Float32"}};
      ptsWriter.addElement("DataArray", attribs);
      Vector<float> vel = ptsWriter.unpackArray(this->u);
      ptsWriter.writeDataArray(vel);
      ptsWriter.closeElement(); // DataArray
    }

    // Point Data 
    ptsWriter.closeElement();
  
    // here's the problem: ParaView's VTK/XML reader is not able to read "raw" bytestreams
    // it parses each character and inevitably sees a > or <, so complains about mismatched
    //   tags, and doesn't seem to point to the right place
    // everyone who complains about raw data in the AppendedData section makes the mistake
    //   of forgetting the 4-byte length - I've got that
    // it's VTK's fault here
    // Jesus, for how inefficient saving 2D particle data is, it's compounded by having to 
    //   do it all in ascii!!! VTK just sucks. That's really what my 6 hours of wasted work
    //   has taught me. It just sucks. Don't use it. Not that anything else is better.
  /*
    printer.OpenElement( "AppendedData" );
    printer.PushAttribute( "encoding", "raw" );
    printer.PushText( " " );
    printer.PushText( "_" );
    uint32_t arry_len = s.size() * sizeof(s[0]);
    char* ptr = (char*)(&arry_len);
    std::cout << "  writing " << arry_len << " bytes to appended data" << std::endl;
    //std::fwrite(&arry_len, sizeof(uint32_t), 1, fp);
    std::fwrite(ptr, sizeof(uint32_t), 1, fp);
    std::fwrite(s.data(), sizeof(s[0]), s.size(), fp);
    printer.PushText( " " );
    printer.CloseElement();	// AppendedData
  */
  
    // Piece 
    ptsWriter.closeElement();
    ptsWriter.finish();
    std::cout << "Wrote " << this->n << " points to " << vtkfn.str() << std::endl;
    return vtkfn.str();
  }

protected:
  // additional state vector
  Vector<S> r;					// thickness/radius
  // in 3D, this is where elong and ug would be

private:
#ifdef USE_GL
  std::shared_ptr<GlState> mgl;
#endif
  float max_strength;
};

