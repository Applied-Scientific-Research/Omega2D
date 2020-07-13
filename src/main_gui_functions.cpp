/*
 * main_gui_functions.cpp - contains helper, non-class related functions used in main_gui.cpp
 *                          to in increase readability
 *
 * (c) 2020 Applied Scientific Reasearch, Inc.
 *          Mark J Stock <markjstock@gmail.com>
 *          Blake B Hillier <blakehillier@mac.com>
 */

#pragma once

#ifdef USE_IMGUI
#include <GLFW/glfw3.h>
#include "imgui/imgui.h"

static void error_callback(int error, const char* description) {
  fprintf(stderr, "Error %d: %s\n", error, description);
}

//static void keyboard_callback(int key, int action) {
//  printf("%d %d\n", key, action);
//}

//
// this is NOT a GLFW callback, but my own, and it needs the
// ImGuiIO data structure for information on the mouse state
//
void mouse_callback(GLFWwindow* /*_thiswin*/,
                    ImGuiIO& io,
                    float*   _cx,
                    float*   _cy,
                    float*   _size) {

  // first, use left-click and drag to move the data
  static bool lbutton_down = false;

  if (io.MouseClicked[0]) lbutton_down = true;
  if (io.MouseReleased[0]) lbutton_down = false;

  if (lbutton_down) {
    //std::cout << "free mouse moved " << io.MouseDelta.x << " " << io.MouseDelta.y << std::endl;
    // do your drag here

    // this worked on Linux:
    //int display_w, display_h;
    //glfwGetFramebufferSize(_thiswin, &display_w, &display_h);
    //(*_cx) -= 2.0f * (*_size) * (float)io.MouseDelta.x / (float)display_w;
    //(*_cy) += 2.0f * (*_size) * (float)io.MouseDelta.y / (float)display_w;

    // this works on a Retina display:
    (*_cx) -= 2.0f * (*_size) * (float)io.MouseDelta.x / io.DisplaySize.x;
    (*_cy) += 2.0f * (*_size) * (float)io.MouseDelta.y / io.DisplaySize.x;
  }

  // then, use scroll wheel to zoom!
  //std::cout << "free mouse wheel " << io.MouseWheel << std::endl;
  if (io.MouseWheel != 0) {
    // do your drag here
    //int display_w, display_h;
    //glfwGetFramebufferSize(_thiswin, &display_w, &display_h);

    // change the size
    const float oldsize = (*_size);
    (*_size) *= std::pow(1.1f, io.MouseWheel);

    // and adjust the center such that the zoom occurs about the pointer!
    //const float ar = (float)display_h / (float)display_w;
    const float ar = io.DisplaySize.y / io.DisplaySize.x;

    // this only scales around world origin
    //(*_cx) += 2.0f * ((float)io.MousePos.x / (float)display_w - 0.5f) * (oldsize - (*_size));
    //(*_cy) += 2.0f * (0.5f - (float)io.MousePos.y / (float)display_h) * (oldsize - (*_size)) * ar;
    (*_cx) += 2.0f * ((float)io.MousePos.x / io.DisplaySize.x - 0.5f) * (oldsize - (*_size));
    (*_cy) += 2.0f * (0.5f - (float)io.MousePos.y / io.DisplaySize.y) * (oldsize - (*_size)) * ar;
  }
}

//
// Helper routine to determine orthographic projection matrix
// given coords at screen center and a measure of size
// Also changes overall pixels-to-length scale
//
void compute_ortho_proj_mat(GLFWwindow*         _thiswin,
                            const float         _cx,
                            const float         _cy,
                            float*              _size,
                            std::vector<float>& _projmat) {

  // track changes in window!
  static int last_w, last_h = -1;

  // get current window size
  int display_w, display_h;
  glfwGetFramebufferSize(_thiswin, &display_w, &display_h);

  // compare window size to previous call
  if (last_h != -1) {
    if (last_h != display_h or last_w != display_w) {
      // window aspect ratio changed, adjust _size
      (*_size) *= sqrt(  ((float)last_h   /(float)last_w   )
                       / ((float)display_h/(float)display_w));
    }
  }

  const float vsx = (*_size);
  const float vsy = (*_size) * (float)display_h / (float)display_w;
  _projmat =
    { 1.0f/vsx, 0.0f,     0.0f, 0.0f,
      0.0f,     1.0f/vsy, 0.0f, 0.0f,
      0.0f,     0.0f,    -1.0f, 0.0f,
     -_cx/vsx, -_cy/vsy,  0.0f, 1.0f };

  // save window size for next call
  last_w = display_w;
  last_h = display_h;
}

//
// resize a window and framebuffer programmatically
//
void resize_to_resolution(GLFWwindow* window, const int new_w, const int new_h) {

  // get framebuffer size
  int fb_w, fb_h;
  glfwGetFramebufferSize(window, &fb_w, &fb_h);
  //std::cout << "Framebuffer size is " << fb_w << " x " << fb_h << std::endl;

  // get window size
  int ws_w, ws_h;
  glfwGetWindowSize(window, &ws_w, &ws_h);
  //std::cout << "Window size is " << ws_w << " x " << ws_h << std::endl;

  // on normal monitors, these numbers should be the same; on retina displays, they may not

  // check and resize
  if (fb_w != new_w or fb_h != new_h) {
    // on a retina display...do anything different?

    glfwSetWindowSize(window, new_w, new_h);
    std::cout << "Resizing window/framebuffer to " << new_w << " x " << new_h << std::endl;
  }
}

void LoadJsonSims(std::vector<nlohmann::json> &sims, std::vector<std::string> &descriptions, const std::string dirPath) {
  std::list<std::string> fileNames = MiniPath::listFiles(dirPath, "*.json");
  std::string sysDelim = MiniPath::getSystemDelim();
  std::cout << "Reading in" << std::endl;
  for(const std::string& s : fileNames) {
    sims.push_back(read_json(dirPath+sysDelim+s));
    descriptions.push_back(sims.back()["description"]);
  }
}

int obj_movement_gui(int &mitem, char* strx, char* stry, char* strrad) {
  // fixed to ground      - this geometry is fixed (attached to inertial)
  // attached to previous - this geometry is attached to the previous geometry
  // according to formula - this geometry is attached to a new moving body
  const char* mitems[] = { "fixed to ground", "attached to previous", "according to formula" };
  int changed = 0;
  static int tmp = -1;
  //const char* mitems[] = { "fixed", "attached to previous", "according to formula", "dynamic" };
  ImGui::Combo("movement", &mitem, mitems, 3);
  if (tmp != mitem) { 
    tmp = mitem;
    changed += 1;
  }
  // show different inputs based on what is selected
  if (mitem == 2) {
    changed += ImGui::InputText("x position", strx, 512);
    ImGui::SameLine();
    ShowHelpMarker("Use C-style expressions, t is time\n+ - / * % ^ ( ) pi e\nabs, sin, cos, tan, exp, log, log10, sqrt, floor, pow");
    changed += ImGui::InputText("y position", stry, 512);
    ImGui::SameLine();
    ShowHelpMarker("Use C-style expressions, t is time\n+ - / * % ^ ( ) pi e\nabs, sin, cos, tan, exp, log, log10, sqrt, floor, pow");
    changed += ImGui::InputText("angular position", strrad, 512);
    ImGui::SameLine();
    ShowHelpMarker("In radians, use C-style expressions, t is time\n+ - / * % ^ ( ) pi e\nabs, sin, cos, tan, exp, log, log10, sqrt, floor, pow");
  }
  
  return changed;
}

// execution starts here
void draw_render_gui(RenderParams &rp) {
  ImGui::ColorEdit3("positive circulation", rp.pos_circ_color);
  ImGui::ColorEdit3("negative circulation", rp.neg_circ_color);
  ImGui::ColorEdit3("feature color",        rp.default_color);
  ImGui::ColorEdit3("background color",     rp.clear_color);
  //ImGui::Checkbox("show origin", &show_origin);
  ImGui::SliderFloat("particle brightness", &(rp.circ_density), 0.01f, 10.0f, "%.2f", 2.0f);
  ImGui::SliderFloat("particle scale", &(rp.vorton_scale), 0.01f, 2.0f, "%.2f", 2.0f);

  if (ImGui::Button("Recenter")) {
    // put everything back to cente
    rp.vcx = -0.5f;
    rp.vcy = 0.0f;
    rp.vsize = 2.0f;
  }
  // add button to recenter on all vorticity?
}

void draw_stats_window(const long int numPanels, const long int numFieldPts, const long int step, const float time,
                       const long int numParticles, bool* showStatsWindow, const int fontSize, const float displayH) {
  //std::cout << "Creating stats window: " << std::endl;
  // there's no way to have this appear in the output png without the rest of the GUI
  const int numrows = 4 + (numPanels > 0 ? 1 : 0) + (numFieldPts > 0 ? 1 : 0);
  // std::cout << "   fontSize: " << fontSize << "\n   numrows: " << numrows << "\n   display_h: " << display_h << std::endl;
#ifdef __APPLE__
  ImGui::SetNextWindowSize(ImVec2(10+fontSize*11,10+1.1*fontSize*numrows));
  ImGui::SetNextWindowPos(ImVec2(10.0f, ((displayH-fontSize*(1.1*numrows))/2)-60.0));
#else
  ImGui::SetNextWindowSize(ImVec2(10+fontSize*11,10+1.1*fontSize*numrows));
  ImGui::SetNextWindowPos(ImVec2(20, displayH-fontSize*(1.1*numrows+1)));
#endif
  ImGuiWindowFlags window_flags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove;
  bool h = ImGui::Begin("Statistics", showStatsWindow, window_flags);
  if (!h) {
    std::cout << "Uh oh" << std::endl;
  }
  ImGui::Text("Step %13ld", step);
  ImGui::Text("Time %13.4f", time);
  if (numPanels > 0) { ImGui::Text("Panels %11ld", numPanels); }
  ImGui::Text("Particles %8ld", numParticles);
  if (numFieldPts > 0) { ImGui::Text("Field Pts %8ld", numFieldPts); }
  ImGui::End();
}

bool draw_welcome_window(const float displayW, const float displayH) {
  bool show = true;
  ImGui::OpenPopup("Welcome!");
  ImGui::SetNextWindowSize(ImVec2(500,300));
  ImGui::SetNextWindowPos(ImVec2(displayW * 0.5f, displayH * 0.5f), ImGuiCond_Always, ImVec2(0.5f,0.5f));
  ImGuiWindowFlags window_flags = ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove;
  if (ImGui::BeginPopupModal("Welcome!", NULL, window_flags)) {
    //ImGui::Begin("Welcome", &show_welcome_window);
    ImGui::TextWrapped("Welcome to Omega2D! Select a simulation from the drop-down, load from a file, or manually set your simulation global properites and add one or more flow, boundary, or measurement structures. Space bar starts and stops the run, Reset clears and loads new simulation properties. Left mouse button drags the frame around, mouse scroll wheel zooms. Save your flow set-up to json or your flow image to png or vtu. Have fun!");
    ImGui::Spacing();
    //if (ImGui::Button("Got it.", ImVec2(120,0))) { show_welcome_window = false; }
    //ImGui::End();
    // const float xwid = ImGui::GetWindowContentRegionWidth();
    if (ImGui::Button("Got it", ImVec2(120,0))) {
      ImGui::CloseCurrentPopup(); 
      show = false;
    }
    ImGui::EndPopup();
  }
  return show;
}
#endif
