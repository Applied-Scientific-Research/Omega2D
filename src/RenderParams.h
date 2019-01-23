/*
 * RenderParams.h - Structure to contain rendering/drawing parameters
 *
 * (c)2019 Applied Scientific Research, Inc.
 *         Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

//
// Class-like struct for all Imgui and OpenGL render parameters
//
struct RenderParams {
  RenderParams() = default;
  ~RenderParams() = default;

  RenderParams(RenderParams const&) = default; //allow copy
  RenderParams(RenderParams&&) = default; //allow move
  RenderParams& operator=(RenderParams const&) = default; //allow copy
  RenderParams& operator=(RenderParams&&) = default; //allow move

  // public-equivalent data

  // window space
  size_t width = 1280;
  size_t height = 720;

  // world space
  float vcx = -0.5f;
  float vcy = 0.0f;
  float vsize = 2.0f;

  // colors
  float pos_circ_color[4] = {207./255.,  47./255.,  47./255., 1.0};
  float neg_circ_color[4] = { 63./255.,  63./255., 255./255., 1.0};
  float default_color[4]  = {204./255., 204./255., 204./255., 1.0};
  float clear_color[4]    = { 15./255.,  15./255.,  15./255., 1.0};

  // other
  float tracer_scale = 0.15;	// non-dimensional
  float tracer_size;		// dimensional
};

