/*
 * RenderParams.h - Structure to contain rendering/drawing parameters
 *
 * (c)2019 Applied Scientific Research, Inc.
 *         Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "imgui/imgui.h"

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
  ImVec4 pos_circ_color = ImColor(207, 47, 47);
  ImVec4 neg_circ_color = ImColor(63, 63, 255);
  ImVec4 default_color  = ImColor(204, 204, 204);
  ImVec4 clear_color    = ImColor(15, 15, 15);

  // other
  float tracer_scale = 0.15;	// non-dimensional
  float tracer_size;		// dimensional
};

