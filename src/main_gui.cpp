/*
 * main_gui.cpp - Driver code for Omega2D + ImGui + Vc vortex particle method
 *                and boundary element method solver, GUI version
 *
 * (c)2017-20 Applied Scientific Research, Inc.
 *            Written by Mark J Stock <markjstock@gmail.com>
 */

#include "FlowFeature.h"
#include "BoundaryFeature.h"
#include "MeasureFeature.h"
#include "Simulation.h"
#include "JsonHelper.h"
#include "Body.h"
#include "RenderParams.h"

#ifdef _WIN32
  // for glad
  #ifndef APIENTRY
    #define APIENTRY __stdcall
  #endif
  // for C++11 stuff that Windows can't get right
  #include <ciso646>
#endif
#include "glad.h"

// header-only immediate-mode GUI
#include "GuiHelper.h"

// header-only png writing
#include "stb/FrameBufferToImage.h"

//#include <GL/gl3w.h>    // This example is using gl3w to access OpenGL
// functions (because it is small). You may use glew/glad/glLoadGen/etc.
// whatever already works for you.
#include <GLFW/glfw3.h>

#include <cstdio>
#include <iostream>
#include <vector>
#include <iomanip>	// for setfill, setw
//#include <fenv.h>

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


// execution starts here

int main(int argc, char const *argv[]) {
  std::cout << std::endl << "Omega2D GUI" << std::endl;

  // Set up vortex particle simulation
  Simulation sim;
  std::vector< std::unique_ptr<FlowFeature> > ffeatures;
  std::vector< std::unique_ptr<BoundaryFeature> > bfeatures;
  std::vector< std::unique_ptr<MeasureFeature> > mfeatures;
  size_t nframes = 0;
  static bool sim_is_running = false;
  static bool begin_single_step = false;

  // placeholder for command-line input file
  std::string command_line_input;

  // Set up primary OpenGL window
  glfwSetErrorCallback(error_callback);
  if (!glfwInit())
    return 1;
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
#if __APPLE__
  const char* glsl_version = "#version 150";
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#else
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
  const char* glsl_version = "#version 130";
#endif
  GLFWwindow* window = glfwCreateWindow(1280, 720, "Omega2D GUI", nullptr, nullptr);
  if (!window) { 
  std::cout << "glfwCreateWindow created " << window << std::endl;
  exit(-1);
  }
  glfwMakeContextCurrent(window);
  glfwSwapInterval(1); // Enable vsync

  //gl3wInit();
  if (!gladLoadGL()) {
    std::cout << "gladLoadGL failed " << std::endl;
    exit(-1);
  }

  // Setup ImGui binding
  ImGui::CreateContext();
  ImGui_ImplGlfw_InitForOpenGL(window, true);
  ImGui_ImplOpenGL3_Init(glsl_version);

  //glfwSetKeyCallback(keyboard_callback);

  //glfwSetWindowCloseCallback(window, window_close_callback);

  // Get and set some IO functions
  ImGuiIO& io = ImGui::GetIO();
  io.IniFilename = ".omega2d.ini";
  std::vector<std::string> recent_json_files;

  // Load Fonts
  // (there is a default font, this is only if you want to change it. see extra_fonts/README.txt for more details)

  //  1. (Optional) Call AddFont*** functions. If you don't call any, the default font will be loaded for you.
  //  2. Call GetTexDataAsAlpha8() or GetTexDataAsRGBA32() to build and retrieve pixels data.
  //  3. Upload the pixels data into a texture within your graphics system.
  //  4. Call SetTexID(my_tex_id); and pass the pointer/identifier to your texture. This value will be passed back to you during rendering to identify the texture.

  // increase the font size
  float fontSize = 20.0f;
  ImFontConfig config;
  //config.PixelSnapH = true;
  config.GlyphOffset = ImVec2(0.0f, 1.0f);
  config.SizePixels = fontSize;
  io.Fonts->AddFontDefault(&config);

  // try using my own font
  //io.Fonts->AddFontFromFileTTF("Roboto-Regular.ttf", 22);
  //io.Fonts->AddFontFromFileTTF("DroidSansMono.ttf", 20);
  //{
    //unsigned char* pixels;
    //int width, height;
    //io.Fonts->GetTexDataAsRGBA32(&pixels, &width, &height);
    //std::cout << "New font rasterized as " << width << " by " << height << std::endl;
  //}

  // a string to hold any error messages
  std::string sim_err_msg;

  // GUI and drawing parameters
  bool export_vtk_this_frame = false;	// write a vtk with the current data
  std::vector<std::string> vtk_out_files; // list of just-output files
  bool draw_this_frame = false;		// draw the frame as soon as its done
  std::string png_out_file;		// the name of the recently-written png
  bool record_all_frames = false;	// save a frame when a new one is ready
  bool show_stats_window = true;
  bool show_welcome_window = true;
  bool show_terminal_window = false;
  bool show_test_window = false;
  bool show_json_input_window = false;
  bool show_file_output_window = false;
  //static bool show_origin = true;
  static bool is_viscous = false;

  // colors and projection matrix for the render view
  RenderParams rparams;
  std::vector<float> gl_projection;
  compute_ortho_proj_mat(window, rparams.vcx, rparams.vcy, &rparams.vsize, gl_projection);

  // adjust some UI settings
  ImGuiStyle& style = ImGui::GetStyle();
  style.Colors[ImGuiCol_WindowBg]              = ImVec4(0.00f, 0.00f, 0.00f, 1.00f);
  style.Colors[ImGuiCol_TitleBg]               = ImVec4(0.27f, 0.27f, 0.54f, 1.00f);
  style.Colors[ImGuiCol_TitleBgActive]         = ImVec4(0.32f, 0.32f, 0.63f, 1.00f);

  // FE_OVERFLOW is triggered in ImGui
  // FE_INVALID is triggered in std::acos
  // FE_DIVBYZERO is triggered somewhere before Reflect.h:512, and only under x86, when switching from particle-only to BEM
  //feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  //feenableexcept(FE_DIVBYZERO);


  // Main loop
  std::cout << "Starting main loop" << std::endl;
  while (!glfwWindowShouldClose(window))
  {
    glfwPollEvents();
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
    
    //
    // Initialize simulation
    //

    if (not sim.is_initialized() and (sim_is_running || begin_single_step)) {

      std::cout << std::endl << "Initializing simulation" << std::endl;

      // initialize particle distributions
      for (auto const& ff: ffeatures) {
        if (ff->is_enabled()) sim.add_particles( ff->init_particles(sim.get_ips()) );
      }

      // initialize solid objects
      for (auto const& bf : bfeatures) {
        if (bf->is_enabled()) sim.add_boundary( bf->get_body(), bf->init_elements(sim.get_ips()) );
      }

      // initialize measurement features
      for (auto const& mf: mfeatures) {
        if (mf->is_enabled()) sim.add_fldpts( mf->init_particles(rparams.tracer_scale*sim.get_ips()), mf->moves() );
      }

      sim.set_initialized();

      // check setup for obvious errors
      sim_err_msg = sim.check_initialization();

      if (sim_err_msg.empty()) {
        // begin the initial calculation at t=0 so we can draw something
        sim.async_first_step();

      } else {
        // the last step had some difficulty
        std::cout << std::endl << "ERROR: " << sim_err_msg;
        // stop the run
        sim_is_running = false;
      }

      begin_single_step = false;
    }

    //
    // Update simulation
    //

    // get results of latest step, if it just completed
    bool is_ready = sim.test_for_new_results();

    // before we start again, write the vtu output
    if (is_ready and export_vtk_this_frame) {
      vtk_out_files = sim.write_vtk();
      export_vtk_this_frame = false;
    }

    // draw a notification box
    if (not vtk_out_files.empty()) {
      static int32_t vtkframect = 0;
      ++vtkframect;

      // draw the notification
      ImGui::SetNextWindowSize(ImVec2(10+fontSize*12, 10+fontSize*(2+vtk_out_files.size())));
      ImGui::SetNextWindowPos(ImVec2(0,0), 0);
      ImGuiWindowFlags window_flags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoResize;
      ImGui::Begin("Vtk written", NULL, window_flags);
      ImGui::Text("Wrote %ld file(s):", vtk_out_files.size());
      for (auto &thisfile : vtk_out_files) {
        ImGui::Text("  %s", thisfile.c_str());
      }
      ImGui::End();

      // make sure this isn't up for too long
      if (vtkframect == 90) {
        vtkframect = 0;
        vtk_out_files.clear();
      }
    }

    // see if we should start a new step
    if (is_ready and (sim_is_running || begin_single_step)) {

      // check flow for blow-up or dynamic errors
      sim_err_msg = sim.check_simulation();

      if (sim_err_msg.empty()) {
        // the last simulation step was fine, OK to continue

        // generate new particles from emitters
        for (auto const& ff : ffeatures) {
          if (ff->is_enabled()) sim.add_particles( ff->step_particles(sim.get_ips()) );
        }
        for (auto const& mf : mfeatures) {
          if (mf->is_enabled()) sim.add_fldpts( mf->step_particles(rparams.tracer_scale*sim.get_ips()), true );
        }

        // begin a new dynamic step: convection and diffusion
        sim.async_step();

      } else {
        // the last step had some difficulty
        std::cout << std::endl << "ERROR: " << sim_err_msg;

        // stop the run
        sim_is_running = false;
      }

      begin_single_step = false;
    }

    // check the error message and display if there is one
    if (not sim_err_msg.empty()) {
      // write a warning/error message
      ImGui::OpenPopup("Simulation error occurred");
      ImGui::SetNextWindowSize(ImVec2(400,200), ImGuiSetCond_FirstUseEver);
      if (ImGui::BeginPopupModal("Simulation error occurred")) {
        ImGui::Spacing();
        ImGui::TextWrapped(sim_err_msg.c_str());
        ImGui::Spacing();
        if (ImGui::Button("Got it.", ImVec2(120,0))) {
          // clear out the error message first
          sim_err_msg.clear();
          ImGui::CloseCurrentPopup();
        }
        ImGui::EndPopup();
      }
    }


    // check mouse for drag and rescaling!
    if (not io.WantCaptureMouse) {
      mouse_callback(window, io, &rparams.vcx, &rparams.vcy, &rparams.vsize);
    }

    // check for keypresses to toggle state
    //if (not io.WantCaptureKeyboard) {
      //keyboard_callback(
    //}

    //
    // The main window
    //
    {

    ImGui::SetNextWindowSize(ImVec2(140+fontSize*24,100+fontSize*12), ImGuiSetCond_FirstUseEver);
    ImGui::SetNextWindowPos(ImVec2(20,20), ImGuiSetCond_FirstUseEver);
    ImGui::Begin("Omega2D");
    ImGui::Spacing();

    // Select pre-populated simulations
    {
      static int sim_item = 0;
      const char* sim_items[] = { "Select a simulation...", "co-rotating vortices", "traveling vortex pair", "asymmetric vortex patch", "flow over circle", "flow over square", "driven cavity" };
      ImGui::Combo("", &sim_item, sim_items, 7);

      float* dt = sim.addr_dt();
      float* fs = sim.addr_fs();
      float* re = sim.addr_re();
      std::shared_ptr<Body> bp;

      switch(sim_item) {
        case 0:
          // nothing
          break;
        case 1:
          // axisymmetrization of a co-rotating vortex pair
          sim.reset();
          sim.clear_bodies();
          bfeatures.clear();
          ffeatures.clear();
          mfeatures.clear();
          *dt = 0.02;
          fs[0] = 0.0; fs[1] = 0.0;
          sim.set_re_for_ips(0.025);
          // generate the vortices
          ffeatures.emplace_back(std::make_unique<VortexBlob>(0.0, 0.5, 1.0, 0.30, 0.1));
          ffeatures.emplace_back(std::make_unique<VortexBlob>(0.0, -0.5, 1.0, 0.30, 0.1));
          is_viscous = false;
          sim.set_diffuse(false);
          rparams.vcx = -0.5;
          rparams.vcy = 0.0;
          rparams.vsize = 2.0;
          rparams.circ_density = 0.35;
          // start it up
          sim_is_running = true;
          // and make sure we don't keep re-entering this
          sim_item = 0;
          break;
        case 2:
          // travelling counter-rotating pair
          sim.reset();
          sim.clear_bodies();
          bfeatures.clear();
          ffeatures.clear();
          mfeatures.clear();
          *dt = 0.02;
          fs[0] = 0.0; fs[1] = 0.0;
          sim.set_re_for_ips(0.025);
          // generate the vortices
          ffeatures.emplace_back(std::make_unique<VortexBlob>(0.0, 0.5, 1.0, 0.30, 0.1));
          ffeatures.emplace_back(std::make_unique<VortexBlob>(0.0, -0.5, -1.0, 0.30, 0.1));
          is_viscous = false;
          sim.set_diffuse(false);
          rparams.vcx = -0.5;
          rparams.vcy = 0.0;
          rparams.vsize = 2.0;
          rparams.circ_density = 0.35;
          // start it up
          sim_is_running = true;
          // and make sure we don't keep re-entering this
          sim_item = 0;
          break;
        case 3:
          // axi-symmetrization of an oval vortex patch
          sim.reset();
          sim.clear_bodies();
          bfeatures.clear();
          ffeatures.clear();
          mfeatures.clear();
          *dt = 0.02;
          fs[0] = 0.0; fs[1] = 0.0;
          sim.set_re_for_ips(0.015);
          // generate the vortices
          ffeatures.emplace_back(std::make_unique<AsymmetricBlob>(0.0, 0.0, 1.0, 1./2., 1./8., 1/32., 90.));
          is_viscous = true;
          sim.set_diffuse(true);
          rparams.vcx = -0.5;
          rparams.vcy = 0.0;
          rparams.vsize = 2.0;
          rparams.circ_density = 0.6;
          // start it up
          sim_is_running = true;
          // and make sure we don't keep re-entering this
          sim_item = 0;
          break;
        case 4:
          // Re=250 flow over a circular cylinder
          sim.reset();
          sim.clear_bodies();
          bfeatures.clear();
          ffeatures.clear();
          mfeatures.clear();
          *dt = 0.02;
          fs[0] = 1.0; fs[1] = 0.0;
          *re = 250.0;
          // generate the boundary
          bp = sim.get_pointer_to_body("ground");
          bfeatures.emplace_back(std::make_unique<SolidCircle>(bp, true, 0.0, 0.0, 1.0));
          is_viscous = true;
          sim.set_diffuse(true);
          rparams.vcx = -0.5;
          rparams.vcy = 0.0;
          rparams.vsize = 2.0;
          rparams.circ_density = 0.35;
          // start it up
          sim_is_running = true;
          // and make sure we don't keep re-entering this
          sim_item = 0;
          break;
        case 5:
          // Re=500 flow over a square
          sim.reset();
          sim.clear_bodies();
          bfeatures.clear();
          ffeatures.clear();
          mfeatures.clear();
          *dt = 0.01;
          fs[0] = 1.0; fs[1] = 0.0;
          *re = 500.0;
          // generate the boundary
          bp = sim.get_pointer_to_body("ground");
          bfeatures.emplace_back(std::make_unique<SolidSquare>(bp, true, 0.0, 0.0, 1.0, 0.0));
          is_viscous = true;
          sim.set_diffuse(true);
          rparams.vcx = -0.5;
          rparams.vcy = 0.0;
          rparams.vsize = 2.0;
          rparams.circ_density = 0.35;
          // start it up
          sim_is_running = true;
          // and make sure we don't keep re-entering this
          sim_item = 0;
          break;
        case 6:
          // driven cavity
          sim.reset();
          sim.clear_bodies();
          bfeatures.clear();
          ffeatures.clear();
          mfeatures.clear();
          *dt = 0.02;
          fs[0] = 0.0; fs[1] = 0.0;
          *re = 250.0;
          // generate the boundary
          bp = sim.get_pointer_to_body("ground");
          bfeatures.emplace_back(std::make_unique<BoundarySegment>(bp, true, 0.0, 0.0, 1.0, 0.0));
          bfeatures.emplace_back(std::make_unique<BoundarySegment>(bp, true, 1.0, 0.0, 1.0, 1.0));
          bfeatures.emplace_back(std::make_unique<BoundarySegment>(bp, true, 1.0, 1.0, 0.0, 1.0, 0.0, -1.0));
          bfeatures.emplace_back(std::make_unique<BoundarySegment>(bp, true, 0.0, 1.0, 0.0, 0.0));
          is_viscous = true;
          sim.set_diffuse(true);
          rparams.vcx = 0.06;
          rparams.vcy = 0.5;
          rparams.vsize = 1.0;
          rparams.circ_density = 0.35;
          // start it up
          sim_is_running = true;
          // and make sure we don't keep re-entering this
          sim_item = 0;
          break;
      } // end switch
    }

    // or load a simulation from a JSON file
    ImGui::SameLine();
    if (ImGui::Button("Load from json", ImVec2(10+fontSize*8,0))) show_json_input_window = true;

    if (show_json_input_window) {
      bool try_it = false;
      static std::string infile = "input.json";

      if (fileIOWindow( try_it, infile, recent_json_files, "Open", {"*.json", "*.*"}, true, ImVec2(200+26*fontSize,300))) {
        show_json_input_window = false;

        if (try_it and !infile.empty()) {
          // remember
          recent_json_files.push_back( infile );

          // stop and clear before loading
          sim.reset();
          bfeatures.clear();
          ffeatures.clear();
          mfeatures.clear();

          // load and report
          read_json(sim, ffeatures, bfeatures, mfeatures, rparams, infile);

          // we have to manually set this variable
          is_viscous = sim.get_diffuse();

          // run one step so we know what we have, or autostart
          if (sim.autostart()) {
            sim_is_running = true;
          } else {
            begin_single_step = true;
          }

          // check and possibly resize the window to match the saved resolution
          resize_to_resolution(window, rparams.width, rparams.height);
        }
      }
    }
    ImGui::Spacing();

    // or load the sim from the command-line (do this once)
    if (argc == 2 and command_line_input.empty()) {

      // stop and clear before loading
      sim.reset();
      bfeatures.clear();
      ffeatures.clear();
      mfeatures.clear();

      command_line_input = argv[1];
      read_json(sim, ffeatures, bfeatures, mfeatures, rparams, command_line_input);

      // we have to manually set this variable
      is_viscous = sim.get_diffuse();

      // run one step so we know what we have, or autostart
      if (sim.autostart()) {
        sim_is_running = true;
      } else {
        begin_single_step = true;
      }

      // check and possibly resize the window to match the saved resolution
      resize_to_resolution(window, rparams.width, rparams.height);

      // we don't need the welcome banner
      show_welcome_window = false;
    }


    //if (ImGui::CollapsingHeader("Simulation globals", ImGuiTreeNodeFlags_DefaultOpen)) {
    if (ImGui::CollapsingHeader("Simulation globals")) {

      // save current versions, so we know which changed
      const float current_re = sim.get_re();
      const float current_dt = sim.get_dt();

      ImGui::Checkbox("Fluid is viscous (diffuses)", &is_viscous);
      ImGui::SameLine();
      ShowHelpMarker("If checked, simulation will add particles to solve diffusion equation; if unchecked, simulation will run fast, but quickly lose accuracy.");

      if (sim.is_initialized() and is_viscous) {
        ImGui::Text("** Until a reset, Time step and Reynolds number move together **");
      }

      // set input widget width for this and the next few
      ImGui::PushItemWidth(-240);

      ImGui::SliderFloat("Time step", sim.addr_dt(), 0.0001f, 0.1f, "%.4f", 2.0f);
      ImGui::SameLine();
      ShowHelpMarker("Adjust how far into the future each step must simulate. Smaller means better accuracy but slower, larger means lower accuracy but faster.");

      if (is_viscous) {
        sim.set_diffuse(true);

        // and let user choose Reynolds number
        ImGui::SliderFloat("Reynolds number", sim.addr_re(), 10.0f, 10000.0f, "%.1f", 2.0f);
        ImGui::SameLine();
        ShowHelpMarker("Reynolds number is the inverse of viscosity; so larger means less viscosity, smaller particles, and longer run time, but more detail.");

        // if Reynolds number or time step change during a run, adjust the other to keep particle spacing constant
        if (sim.is_initialized()) {
          if (current_re != sim.get_re()) {
            // change dt
            *(sim.addr_dt()) = current_dt * sim.get_re() / current_re;
          } else if (current_dt != sim.get_dt()) {
            // change Re
            *(sim.addr_re()) = current_re * sim.get_dt() / current_dt;
          }
        }

        ImGui::Text("Particle spacing for these settings is %g", sim.get_ips());

      } else {
        static float my_ips = 0.03141;
        ImGui::SliderFloat("Particle spacing", &my_ips, 0.001f, 0.1f, "%.3f", 2.0f);
        ImGui::SameLine();
        ShowHelpMarker("Sets the average size and distance between particles.");
        // change underlying Re when this changes
        sim.set_re_for_ips(my_ips);
        my_ips = sim.get_ips();
      }
      ImGui::InputFloat2("Flow vector", sim.addr_fs());
      ImGui::SameLine();
      ShowHelpMarker("Also called freestream, this is a uniform wind blowing along this vector.");

      ImGui::PopItemWidth();

      // set stop/pause conditions
      bool use_step_pause = sim.using_max_steps();
      ImGui::Checkbox("Pause at step", &use_step_pause);
      if (use_step_pause) {
        int pause_step = sim.get_max_steps();
        ImGui::SameLine();
        ImGui::InputInt(" ", &pause_step);
        sim.set_max_steps((size_t)pause_step);
      } else {
        sim.unset_max_steps();
      }
      bool use_time_pause = sim.using_end_time();
      ImGui::Checkbox("Pause at time", &use_time_pause);
      if (use_time_pause) {
        float pause_time = sim.get_end_time();
        ImGui::SameLine();
        ImGui::InputFloat(" ", &pause_time);
        //ImGui::InputFloat("Pause at time", float* v, float step = 0.0f, float step_fast = 0.0f, int decimal_precision = -1, ImGuiInputTextFlags extra_flags = 0);
        sim.set_end_time((double)pause_time);
      } else {
        sim.unset_end_time();
      }
    }

    ImGui::Spacing();
    //if (ImGui::CollapsingHeader("Flow structures", ImGuiTreeNodeFlags_DefaultOpen)) {
    if (ImGui::CollapsingHeader("Startup structures")) {

      if (ffeatures.size() + bfeatures.size() == 0) {
        ImGui::Text("Add flow or boundry features (like vortex blobs and solid objects) here, then click RUN.");
      }


      ImGui::Spacing();

      // button and modal window for adding new flow structures
      if (ImGui::Button("Add flow")) ImGui::OpenPopup("New flow structure");
      ImGui::SetNextWindowSize(ImVec2(400,200), ImGuiSetCond_FirstUseEver);
      if (ImGui::BeginPopupModal("New flow structure"))
      {
        static int item = 1;
        const char* items[] = { "single particle", "round vortex blob", "asymmetric vortex blob", "block of vorticity", "random particles", "particle emitter" };
        ImGui::Combo("type", &item, items, 6);

        static float xc[2] = {0.0f, 0.0f};
        static float rad = 5.0 * sim.get_ips();
        static float soft = sim.get_ips();
        static float str = 1.0f;
        static int npart = 100;
        static float xs[2] = {2.0f, 2.0f};
        static float strlo = -1.0f;
        static float strhi = 1.0f;

        // always ask for center
        ImGui::InputFloat2("center", xc);

        // show different inputs based on what is selected
        switch(item) {
          case 0: {
            // a single vortex particle
            ImGui::SliderFloat("strength", &str, -1.0f, 1.0f, "%.4f");
            ImGui::TextWrapped("This feature will add 1 particle");
            if (ImGui::Button("Add single particle")) {
              ffeatures.emplace_back(std::make_unique<SingleParticle>(xc[0], xc[1], str));
              std::cout << "Added " << (*ffeatures.back()) << std::endl;
              ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            // it would be nice to be able to put this all in
            //SingleParticle::draw_creation_gui();
            } break;

         case 1: {
            // a blob of multiple vorticies
            ImGui::SliderFloat("strength", &str, -5.0f, 5.0f, "%.4f");
            ImGui::SliderFloat("radius", &rad, sim.get_ips(), 1.0f, "%.4f");
            ImGui::SliderFloat("softness", &soft, sim.get_ips(), 1.0f, "%.4f");
            ImGui::TextWrapped("This feature will add about %d particles", (int)(0.785398175*std::pow((2 * rad + soft) / sim.get_ips(), 2)));
            if (ImGui::Button("Add vortex blob")) {
              ffeatures.emplace_back(std::make_unique<VortexBlob>(xc[0], xc[1], str, rad, soft));
              std::cout << "Added " << (*ffeatures.back()) << std::endl;
              ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            } break;

          case 2: {
            // an asymmetric blob of multiple vorticies
            static float minrad = 2.5 * sim.get_ips();
            static float rotdeg = 90.0f;
            ImGui::SliderFloat("strength", &str, -5.0f, 5.0f, "%.4f");
            ImGui::SliderFloat("major radius", &rad, sim.get_ips(), 1.0f, "%.4f");
            ImGui::SliderFloat("minor radius", &minrad, sim.get_ips(), 1.0f, "%.4f");
            ImGui::SliderFloat("softness", &soft, sim.get_ips(), 1.0f, "%.4f");
            ImGui::SliderFloat("orientation", &rotdeg, 0.0f, 179.0f, "%.0f");
            ImGui::TextWrapped("This feature will add about %d particles", (int)(0.785398175*std::pow((2*rad+soft)/sim.get_ips(), 2)));
            if (ImGui::Button("Add asymmetric vortex blob")) {
              ffeatures.emplace_back(std::make_unique<AsymmetricBlob>(xc[0], xc[1], str, rad, minrad, soft, rotdeg));
              std::cout << "Added " << (*ffeatures.back()) << std::endl;
              ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            } break;

          case 3: {
            // random particles in a rectangle
            ImGui::SliderFloat("strength", &str, -5.0f, 5.0f, "%.4f");
            ImGui::SliderFloat2("box size", xs, 0.01f, 10.0f, "%.4f", 2.0f);
            ImGui::TextWrapped("This feature will add %d particles", (int)(xs[0]*xs[1]/std::pow(sim.get_ips(),2)));
            if (ImGui::Button("Add block of vorticies")) {
              ffeatures.emplace_back(std::make_unique<UniformBlock>(xc[0], xc[1], xs[0], xs[1], str));
              std::cout << "Added " << (*ffeatures.back()) << std::endl;
              ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            } break;

          case 4: {
            // random particles in a rectangle
            ImGui::SliderInt("number", &npart, 1, 10000);
            ImGui::SliderFloat2("box size", xs, 0.01f, 10.0f, "%.4f", 2.0f);
            ImGui::DragFloatRange2("strength range", &strlo, &strhi, 0.001f, -0.1f, 0.1f);
            ImGui::TextWrapped("This feature will add %d particles", npart);
            if (ImGui::Button("Add random vorticies")) {
              ffeatures.emplace_back(std::make_unique<BlockOfRandom>(xc[0], xc[1], xs[0], xs[1], strlo, strhi, npart));
              std::cout << "Added " << (*ffeatures.back()) << std::endl;
              ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            } break;

          case 5: {
            // create a particle emitter
            static float estr = 0.1f;
            ImGui::SliderFloat("strength", &estr, -0.1f, 0.1f, "%.4f");
            ImGui::TextWrapped("This feature will add 1 particle per time step");
            if (ImGui::Button("Add particle emitter")) {
              // this is C++11
              ffeatures.emplace_back(std::make_unique<ParticleEmitter>(xc[0], xc[1], estr));
              std::cout << "Added " << (*ffeatures.back()) << std::endl;
              ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            } break;
        }

        if (ImGui::Button("Cancel", ImVec2(120,0))) { ImGui::CloseCurrentPopup(); }
        ImGui::EndPopup();
      }


      // button and modal window for adding new boundary objects
      ImGui::SameLine();
      if (ImGui::Button("Add boundary")) ImGui::OpenPopup("New boundary structure");
      ImGui::SetNextWindowSize(ImVec2(400,275), ImGuiSetCond_FirstUseEver);
      if (ImGui::BeginPopupModal("New boundary structure"))
      {
        // define movement first
        static int mitem = 0;
        const char* mitems[] = { "fixed to ground", "attached to previous", "according to formula" };
        //const char* mitems[] = { "fixed", "attached to previous", "according to formula", "dynamic" };
        ImGui::Combo("movement", &mitem, mitems, 3);
        static char strx[512] = "0.0*t";
        static char stry[512] = "0.0*t";
        static char strrad[512] = "0.0*t";

        // show different inputs based on what is selected
        switch(mitem) {
          case 0:
            // this geometry is fixed (attached to inertial)
            break;
          case 1:
            // this geometry is attached to the previous geometry
            break;
          case 2:
            // this geometry is attached to a new moving body
            ImGui::InputText("x position", strx, 512);
            ImGui::SameLine();
            ShowHelpMarker("Use C-style expressions, t is time\n+ - / * % ^ ( ) pi e\nabs, sin, cos, tan, exp, log, log10, sqrt, floor, pow");
            ImGui::InputText("y position", stry, 512);
            ImGui::SameLine();
            ShowHelpMarker("Use C-style expressions, t is time\n+ - / * % ^ ( ) pi e\nabs, sin, cos, tan, exp, log, log10, sqrt, floor, pow");
            ImGui::InputText("angular position", strrad, 512);
            ImGui::SameLine();
            ShowHelpMarker("In radians, use C-style expressions, t is time\n+ - / * % ^ ( ) pi e\nabs, sin, cos, tan, exp, log, log10, sqrt, floor, pow");
            break;
        }

        // define geometry second
        static int item = 0;
        const char* items[] = { "circle", "square", "oval", "rectangle", "segment" };
        ImGui::Spacing();
        ImGui::Combo("geometry type", &item, items, 5);

        static bool external_flow = true;

        static float xc[2] = {0.0f, 0.0f};
        static float rotdeg = 0.0f;
        static float circdiam = 1.0;
        static float sqside = 1.0;


        // show different inputs based on what is selected
        switch(item) {
          case 0: {
            // create a circular boundary
            ImGui::Checkbox("Object is in flow", &external_flow);
            ImGui::SameLine();
            ShowHelpMarker("Keep checked if object is immersed in flow,\nuncheck if flow is inside of object");
            ImGui::InputFloat2("center", xc);
            ImGui::SliderFloat("diameter", &circdiam, 0.01f, 10.0f, "%.4f", 2.0);
            ImGui::TextWrapped("This feature will add a solid circular boundary centered at the given coordinates");
            if (ImGui::Button("Add circular boundary")) {
              std::shared_ptr<Body> bp;
              switch(mitem) {
                case 0:
                  // this geometry is fixed (attached to inertial)
                  bp = sim.get_pointer_to_body("ground");
                  break;
                case 1:
                  // this geometry is attached to the previous geometry (or ground)
                  bp = sim.get_last_body();
                  break;
                case 2:
                  // this geometry is attached to a new moving body
                  bp = std::make_shared<Body>();
                  bp->set_pos(0, std::string(strx));
                  bp->set_pos(1, std::string(stry));
                  bp->set_rot(std::string(strrad));
                  bp->set_name("circular cylinder");
                  sim.add_body(bp);
                  break;
              }
              bfeatures.emplace_back(std::make_unique<SolidCircle>(bp, external_flow, xc[0], xc[1], circdiam));
              std::cout << "Added " << (*bfeatures.back()) << std::endl;
              ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            } break;

          case 1: {
            // create a square/rectangle boundary
            ImGui::Checkbox("Object is in flow", &external_flow);
            ImGui::SameLine();
            ShowHelpMarker("Keep checked if object is immersed in flow,\nuncheck if flow is inside of object");
            ImGui::InputFloat2("center", xc);
            ImGui::SliderFloat("side length", &sqside, 0.1f, 10.0f, "%.4f");
            ImGui::SliderFloat("orientation", &rotdeg, 0.0f, 89.0f, "%.0f");
            //ImGui::SliderAngle("orientation", &rotdeg);
            ImGui::TextWrapped("This feature will add a solid square boundary centered at the given coordinates");
            if (ImGui::Button("Add square boundary")) {
              std::shared_ptr<Body> bp;
              switch(mitem) {
                case 0:
                  // this geometry is fixed (attached to inertial)
                  bp = sim.get_pointer_to_body("ground");
                  break;
                case 1:
                  // this geometry is attached to the previous geometry (or ground)
                  bp = sim.get_last_body();
                  break;
                case 2:
                  // this geometry is attached to a new moving body
                  bp = std::make_shared<Body>();
                  bp->set_pos(0, std::string(strx));
                  bp->set_pos(1, std::string(stry));
                  bp->set_rot(std::string(strrad));
                  bp->set_name("square cylinder");
                  sim.add_body(bp);
                  break;
              }
              bfeatures.emplace_back(std::make_unique<SolidSquare>(bp, external_flow, xc[0], xc[1], sqside, rotdeg));
              std::cout << "Added " << (*bfeatures.back()) << std::endl;
              ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            } break;

          case 2: {
            // create an oval boundary
            ImGui::Checkbox("Object is in flow", &external_flow);
            ImGui::SameLine();
            ShowHelpMarker("Keep checked if object is immersed in flow,\nuncheck if flow is inside of object");
            ImGui::InputFloat2("center", xc);
            static float minordiam = 0.5;
            ImGui::SliderFloat("major diameter", &circdiam, 0.01f, 10.0f, "%.4f", 2.0);
            ImGui::SliderFloat("minor diameter", &minordiam, 0.01f, 10.0f, "%.4f", 2.0);
            ImGui::SliderFloat("orientation", &rotdeg, 0.0f, 179.0f, "%.0f");
            ImGui::TextWrapped("This feature will add a solid oval boundary centered at the given coordinates");
            if (ImGui::Button("Add oval boundary")) {
              std::shared_ptr<Body> bp;
              switch(mitem) {
                case 0:
                  // this geometry is fixed (attached to inertial)
                  bp = sim.get_pointer_to_body("ground");
                  break;
                case 1:
                  // this geometry is attached to the previous geometry (or ground)
                  bp = sim.get_last_body();
                  break;
                case 2:
                  // this geometry is attached to a new moving body
                  bp = std::make_shared<Body>();
                  bp->set_pos(0, std::string(strx));
                  bp->set_pos(1, std::string(stry));
                  bp->set_rot(std::string(strrad));
                  bp->set_name("oval cylinder");
                  sim.add_body(bp);
                  break;
              }
              bfeatures.emplace_back(std::make_unique<SolidOval>(bp, external_flow, xc[0], xc[1], circdiam, minordiam, rotdeg));
              std::cout << "Added " << (*bfeatures.back()) << std::endl;
              ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            } break;

          case 3: {
            // create a rectangle boundary
            ImGui::Checkbox("Object is in flow", &external_flow);
            ImGui::SameLine();
            ShowHelpMarker("Keep checked if object is immersed in flow,\nuncheck if flow is inside of object");
            ImGui::InputFloat2("center", xc);
            static float rectside = 0.5;
            ImGui::SliderFloat("horizontal size", &sqside, 0.1f, 10.0f, "%.4f");
            ImGui::SliderFloat("vertical size", &rectside, 0.1f, 10.0f, "%.4f");
            ImGui::SliderFloat("orientation", &rotdeg, 0.0f, 89.0f, "%.0f");
            //ImGui::SliderAngle("orientation", &rotdeg);
            ImGui::TextWrapped("This feature will add a solid rectangular boundary centered at the given coordinates");
            if (ImGui::Button("Add rectangular boundary")) {
              std::shared_ptr<Body> bp;
              switch(mitem) {
                case 0:
                  // this geometry is fixed (attached to inertial)
                  bp = sim.get_pointer_to_body("ground");
                  break;
                case 1:
                  // this geometry is attached to the previous geometry (or ground)
                  bp = sim.get_last_body();
                  break;
                case 2:
                  // this geometry is attached to a new moving body
                  bp = std::make_shared<Body>();
                  bp->set_pos(0, std::string(strx));
                  bp->set_pos(1, std::string(stry));
                  bp->set_rot(std::string(strrad));
                  bp->set_name("rectangular cylinder");
                  sim.add_body(bp);
                  break;
              }
              bfeatures.emplace_back(std::make_unique<SolidRect>(bp, external_flow, xc[0], xc[1], sqside, rectside, rotdeg));
              std::cout << "Added " << (*bfeatures.back()) << std::endl;
              ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            } break;

          case 4: {
            // create a straight boundary segment
            static float xe[2] = {1.0f, 0.0f};
            static float tangbc = 0.0;
            ImGui::InputFloat2("start", xc);
            ImGui::InputFloat2("end", xe);
            ImGui::SliderFloat("force tangential flow", &tangbc, -2.0f, 2.0f, "%.1f");
            ImGui::TextWrapped("This feature will add a solid boundary segment from start to end, where fluid is on the left when marching from start to end, and positive tangential flow is as if segment is moving along vector from start to end. Make sure enough segments are created to fully enclose a volume.");
            if (ImGui::Button("Add boundary segment")) {
              std::shared_ptr<Body> bp;
              switch(mitem) {
                case 0:
                  // this geometry is fixed (attached to inertial)
                  bp = sim.get_pointer_to_body("ground");
                  break;
                case 1:
                  // this geometry is attached to the previous geometry (or ground)
                  bp = sim.get_last_body();
                  break;
                case 2:
                  // this geometry is attached to a new moving body
                  bp = std::make_shared<Body>();
                  bp->set_pos(0, std::string(strx));
                  bp->set_pos(1, std::string(stry));
                  bp->set_rot(std::string(strrad));
                  bp->set_name("segmented boundary");
                  sim.add_body(bp);
                  break;
              }
              bfeatures.emplace_back(std::make_unique<BoundarySegment>(bp, true, xc[0], xc[1], xe[0], xe[1], 0.0, tangbc));
              std::cout << "Added " << (*bfeatures.back()) << std::endl;
              ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            } break;
        }

        if (ImGui::Button("Cancel", ImVec2(120,0))) { ImGui::CloseCurrentPopup(); }
        ImGui::EndPopup();
      }


      // button and modal window for adding new measurement objects
      ImGui::SameLine();
      if (ImGui::Button("Add measurement")) ImGui::OpenPopup("New measurement structure");
      ImGui::SetNextWindowSize(ImVec2(400,200), ImGuiSetCond_FirstUseEver);
      if (ImGui::BeginPopupModal("New measurement structure"))
      {
        static int item = 0;
        const char* items[] = { "single point/tracer", "streakline", "circle of tracers", "line of tracers", "measurement line", "measurement grid" };
        ImGui::Combo("type", &item, items, 6);

        static float xc[2] = {0.0f, 0.0f};
        static float xf[2] = {0.0f, 1.0f};
        static bool is_lagrangian = true;
        static float rad = 5.0 * sim.get_ips();

        // show different inputs based on what is selected
        switch(item) {
          case 0: {
            // a single measurement point
            ImGui::InputFloat2("position", xc);
            ImGui::Checkbox("Point follows flow", &is_lagrangian);
            ImGui::TextWrapped("This feature will add 1 point");
            if (ImGui::Button("Add single point")) {
              mfeatures.emplace_back(std::make_unique<SinglePoint>(xc[0], xc[1], is_lagrangian));
              std::cout << "Added " << (*mfeatures.back()) << std::endl;
              ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            } break;
          case 1: {
            // a tracer emitter
            ImGui::InputFloat2("position", xc);
            ImGui::TextWrapped("This feature will add 1 tracer emitter");
            if (ImGui::Button("Add streakline")) {
              mfeatures.emplace_back(std::make_unique<TracerEmitter>(xc[0], xc[1]));
              std::cout << "Added " << (*mfeatures.back()) << std::endl;
              ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            } break;
          case 2: {
            // a tracer circle
            ImGui::InputFloat2("center", xc);
            ImGui::SliderFloat("radius", &rad, sim.get_ips(), 1.0f, "%.4f");
            ImGui::TextWrapped("This feature will add about %d field points",
                               (int)(0.785398175*std::pow(2*rad/(rparams.tracer_scale*sim.get_ips()), 2)));
            if (ImGui::Button("Add circle of tracers")) {
              mfeatures.emplace_back(std::make_unique<TracerBlob>(xc[0], xc[1], rad));
              std::cout << "Added " << (*mfeatures.back()) << std::endl;
              ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            } break;
          case 3: {
            // a tracer line
            ImGui::InputFloat2("start", xc);
            ImGui::InputFloat2("finish", xf);
            ImGui::TextWrapped("This feature will add about %d field points",
                               1+(int)(std::sqrt(std::pow(xf[0]-xc[0],2)+std::pow(xf[1]-xc[1],2))/(rparams.tracer_scale*sim.get_ips())));
            if (ImGui::Button("Add line of tracers")) {
              mfeatures.emplace_back(std::make_unique<TracerLine>(xc[0], xc[1], xf[0], xf[1]));
              std::cout << "Added " << (*mfeatures.back()) << std::endl;
              ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            } break;
          case 4: {
            // a static, measurement line
            ImGui::InputFloat2("start", xc);
            ImGui::InputFloat2("finish", xf);
            ImGui::TextWrapped("This feature will add about %d field points",
                               1+(int)(std::sqrt(std::pow(xf[0]-xc[0],2)+std::pow(xf[1]-xc[1],2))/(rparams.tracer_scale*sim.get_ips())));
            if (ImGui::Button("Add line of measurement points")) {
              mfeatures.emplace_back(std::make_unique<MeasurementLine>(xc[0], xc[1], xf[0], xf[1]));
              std::cout << "Added " << (*mfeatures.back()) << std::endl;
              ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            } break;
          case 5: {
            // a static grid of measurement points
            ImGui::InputFloat2("start", xc);
            ImGui::InputFloat2("finish", xf);
            ImGui::SliderFloat("dx", &rad, sim.get_ips(), 1.0f, "%.4f");
            ImGui::TextWrapped("This feature will add about %d field points",
                               1+(int)((xf[0]-xc[0])*(xf[1]-xc[1])/(rad*rad)));
            if (ImGui::Button("Add grid of measurement points")) {
              mfeatures.emplace_back(std::make_unique<GridPoints>(xc[0], xc[1], xf[0], xf[1], rad));
              std::cout << "Added " << (*mfeatures.back()) << std::endl;
              ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            } break;
        }

        if (ImGui::Button("Cancel", ImVec2(120,0))) { ImGui::CloseCurrentPopup(); }
        ImGui::EndPopup();
      }


      ImGui::Spacing();
      int buttonIDs = 10;

      // list existing flow features here
      int del_this_item = -1;
      for (int i=0; i<(int)ffeatures.size(); ++i) {

        ImGui::PushID(++buttonIDs);
        ImGui::Checkbox("", ffeatures[i]->addr_enabled());
        ImGui::PopID();
        if (ffeatures[i]->is_enabled()) {
          ImGui::SameLine();
          ImGui::Text("%s", ffeatures[i]->to_string().c_str());
        } else {
          ImGui::SameLine();
          ImGui::TextColored(ImVec4(0.5f,0.5f,0.5f,1.0f), "%s", ffeatures[i]->to_string().c_str());
        }

        // add a "remove" button at the end of the line (so it's not easy to accidentally hit)
        ImGui::SameLine();
        ImGui::PushID(++buttonIDs);
        if (ImGui::SmallButton("remove")) del_this_item = i;
        ImGui::PopID();

        //ImGui::SameLine();
        //ImGui::PushID(++buttonIDs);
        //if (ImGui::SmallButton("edit", ImVec2(60,0))) edit_this_item = i;
        //ImGui::PopID();
      }
      if (del_this_item > -1) {
        std::cout << "Asked to delete flow feature " << del_this_item << std::endl;
        ffeatures.erase(ffeatures.begin()+del_this_item);
      }

      // list existing boundary features here
      int del_this_bdry = -1;
      for (int i=0; i<(int)bfeatures.size(); ++i) {

        ImGui::PushID(++buttonIDs);
        ImGui::Checkbox("", bfeatures[i]->addr_enabled());
        ImGui::PopID();
        if (bfeatures[i]->is_enabled()) {
          ImGui::SameLine();
          ImGui::Text("%s", bfeatures[i]->to_string().c_str());
        } else {
          ImGui::SameLine();
          ImGui::TextColored(ImVec4(0.5f,0.5f,0.5f,1.0f), "%s", bfeatures[i]->to_string().c_str());
        }

        // add a "remove" button at the end of the line (so it's not easy to accidentally hit)
        ImGui::SameLine();
        ImGui::PushID(++buttonIDs);
        if (ImGui::SmallButton("remove")) del_this_bdry = i;
        ImGui::PopID();
      }
      if (del_this_bdry > -1) {
        std::cout << "Asked to delete boundary feature " << del_this_bdry << std::endl;
        bfeatures.erase(bfeatures.begin()+del_this_bdry);
      }

      // list existing measurement features here
      int del_this_measure = -1;
      for (int i=0; i<(int)mfeatures.size(); ++i) {

        ImGui::PushID(++buttonIDs);
        ImGui::Checkbox("", mfeatures[i]->addr_enabled());
        ImGui::PopID();
        if (mfeatures[i]->is_enabled()) {
          ImGui::SameLine();
          ImGui::Text("%s", mfeatures[i]->to_string().c_str());
        } else {
          ImGui::SameLine();
          ImGui::TextColored(ImVec4(0.5f,0.5f,0.5f,1.0f), "%s", mfeatures[i]->to_string().c_str());
        }

        // add a "remove" button at the end of the line (so it's not easy to accidentally hit)
        ImGui::SameLine();
        ImGui::PushID(++buttonIDs);
        if (ImGui::SmallButton("remove")) del_this_measure = i;
        ImGui::PopID();
      }
      if (del_this_measure > -1) {
        std::cout << "Asked to delete measurement feature " << del_this_measure << std::endl;
        mfeatures.erase(mfeatures.begin()+del_this_measure);
      }

      //if (ffeatures.size() + bfeatures.size() + mfeatures.size() == 0) {
      //  ImGui::Text("none");
      //}

    } // end structure entry


    // Rendering parameters, under a header
    ImGui::Spacing();
    if (ImGui::CollapsingHeader("Rendering controls")) {
      ImGui::ColorEdit3("positive circulation", rparams.pos_circ_color);
      ImGui::ColorEdit3("negative circulation", rparams.neg_circ_color);
      ImGui::ColorEdit3("feature color",        rparams.default_color);
      ImGui::ColorEdit3("background color",     rparams.clear_color);
      //ImGui::Checkbox("show origin", &show_origin);
      ImGui::SliderFloat("particle brightness", &(rparams.circ_density), 0.01f, 10.0f, "%.2f", 2.0f);
      ImGui::SliderFloat("particle scale", &(rparams.vorton_scale), 0.01f, 2.0f, "%.2f", 2.0f);

      if (ImGui::Button("Recenter")) {
        // put everything back to center
        rparams.vcx = -0.5f;
        rparams.vcy = 0.0f;
        rparams.vsize = 2.0f;
      }

      // add button to recenter on all vorticity?
    }

    // Solver parameters, under its own header
    ImGui::Spacing();
    if (ImGui::CollapsingHeader("Solver parameters (advanced)")) {
      sim.draw_advanced();
    }

    // Output buttons, under a header
    ImGui::Spacing();
    if (ImGui::CollapsingHeader("Save output")) {

      // save the simulation to a JSON or VTK file
      ImGui::Spacing();
      if (ImGui::Button("Save setup to JSON", ImVec2(20+12*fontSize,0))) show_file_output_window = true;
      ImGui::SameLine();
      // PNG output of the render frame
      if (ImGui::Button("Save screenshot to PNG", ImVec2(20+12*fontSize,0))) draw_this_frame = true;

      // next line: VTK output and record
      if (ImGui::Button("Save parts to VTU", ImVec2(20+12*fontSize,0))) export_vtk_this_frame = true;
      ImGui::SameLine();
      if (record_all_frames) {
        if (ImGui::Button("STOP RECORDING", ImVec2(20+12*fontSize,0))) {
          record_all_frames = false;
          sim_is_running = false;
        }
      } else {
        if (ImGui::Button("RECORD to PNG", ImVec2(20+12*fontSize,0))) {
          record_all_frames = true;
          sim_is_running = true;
        }
      }
    }

    if (show_file_output_window) {
      bool try_it = false;
      static std::string outfile = "output.json";

      if (fileIOWindow( try_it, outfile, recent_json_files, "Save", {"*.json", "*.*"}, false, ImVec2(200+26*fontSize,300))) {
        show_file_output_window = false;

        if (try_it) {
          // remember
          recent_json_files.push_back( outfile );

          // retrieve window sizes
          glfwGetWindowSize(window, &rparams.width, &rparams.height);

          // write and echo
          write_json(sim, ffeatures, bfeatures, mfeatures, rparams, outfile);
          std::cout << std::endl << "Wrote simulation to " << outfile << std::endl;
        }
      }
    }



    nframes++;

    // check vs. end conditions, if present
    if (sim.test_vs_stop_async()) {
      // just pause sim
      sim_is_running = false;

      // fully quit if asked
      if (sim.quitonstop()) {
        if (is_ready) {
          // this simulation step has finished, write png and exit
          draw_this_frame = true;

          // tell glfw to close the window next time around
          glfwSetWindowShouldClose(window, GLFW_TRUE);
        }
      }
    }

    // all the other stuff
    {
      ImGui::Spacing();
      ImGui::Separator();
      ImGui::Spacing();

      if (sim_is_running) {
        //ImGui::Text("Simulation is running...step = %ld, time = %g", sim.get_nstep(), sim.get_time());
        if (ImGui::Button("PAUSE", ImVec2(200,20+fontSize))) sim_is_running = false;
        // space bar pauses
        if (ImGui::IsKeyPressed(32) and not show_file_output_window) sim_is_running = false;
      } else {
        //ImGui::Text("Simulation is not running, step = %ld, time = %g", sim.get_nstep(), sim.get_time());
        if (ImGui::Button("RUN", ImVec2(200,20+fontSize))) sim_is_running = true;
        ImGui::SameLine();
        if (ImGui::Button("Step", ImVec2(120,0))) begin_single_step = true;
        // and space bar resumes
        if (ImGui::IsKeyPressed(32) and not show_file_output_window) sim_is_running = true;
      }
      ImGui::SameLine();
      if (ImGui::Button("Reset", ImVec2(120,0))) {
        std::cout << std::endl << "Reset requested" << std::endl;
        // remove all particles and reset timer
        sim.reset();
        std::cout << "Reset complete" << std::endl;
      }

      //ImGui::Spacing();
      //ImGui::Separator();
      /*
      ImGui::Text("Open additional windows");
      if (ImGui::Button("Plot statistics")) show_stats_window ^= 1;
      ImGui::SameLine();
      if (ImGui::Button("Show terminal output")) show_terminal_window ^= 1;
      ImGui::SameLine();
      */

      if (ImGui::Button("ImGui Samples")) show_test_window ^= 1;
      // use ASCII table for number: http://www.asciitable.com/
      // but use CAPITAL letter for a letter, jesus, really?!?
      if (ImGui::IsKeyPressed(84) and not show_file_output_window) show_test_window ^= 1;

      //ImGui::Text("Draw frame rate: %.2f ms/frame (%.1f FPS)",
      //            1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);

      //ImGui::Text("Number of panels: %ld  particles: %ld  field points: %ld",
      //            sim.get_npanels(), sim.get_nparts(), sim.get_nfldpts());
    }


    // done drawing the UI window
    ImGui::End();
    }

    // Rendering
    int display_w, display_h;
    glfwGetFramebufferSize(window, &display_w, &display_h);
    glViewport(0, 0, display_w, display_h);
    glClearColor(rparams.clear_color[0], rparams.clear_color[1], rparams.clear_color[2], rparams.clear_color[3]);
    glClear(GL_COLOR_BUFFER_BIT);

    // Show the welcome window
    if (show_welcome_window) {
      // std::cout << "Welocome!" << std::endl;
      ImGui::OpenPopup("Welcome!");
      ImGui::SetNextWindowSize(ImVec2(500,300));
      ImGui::SetNextWindowPos(ImVec2(io.DisplaySize.x * 0.5f, io.DisplaySize.y * 0.5f), ImGuiCond_Always, ImVec2(0.5f,0.5f));
      ImGuiWindowFlags window_flags = ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove;
      if (ImGui::BeginPopupModal("Welcome!", NULL, window_flags)) {
        //ImGui::Begin("Welcome", &show_welcome_window);
        ImGui::TextWrapped("Welcome to Omega2D! Select a simulation from the drop-down, load from a file, or manually set your simulation global properites and add one or more flow, boundary, or measurement structures. Space bar starts and stops the run, Reset clears and loads new simulation properties. Left mouse button drags the frame around, mouse scroll wheel zooms. Save your flow set-up to json or your flow image to png or vtu. Have fun!");
        ImGui::Spacing();
        //if (ImGui::Button("Got it.", ImVec2(120,0))) { show_welcome_window = false; }
        //ImGui::End();
        // const float xwid = ImGui::GetWindowContentRegionWidth();
        if (ImGui::Button("Got it", ImVec2(120,0))) { ImGui::CloseCurrentPopup(); show_welcome_window = false; }
        ImGui::EndPopup();
      }
    }

    // Show the simulation stats in the corner
    //if (nframes < 10){ std::cout << "show_stats_window: " << show_stats_window << std::endl; }
    if (show_stats_window) {
      // std::cout << "Creating stats window: " << std::endl;
      // there's no way to have this appear in the output png without the rest of the GUI
      const int numrows = 4 + (sim.get_npanels()>0 ? 1 : 0) + (sim.get_nfldpts()>0 ? 1 : 0);
      // std::cout << "   fontSize: " << fontSize << "\n   numrows: " << numrows << "\n   display_h: " << display_h << std::endl;
      ImGui::SetNextWindowSize(ImVec2(10+fontSize*11, 10+1.1*fontSize*numrows));
      ImGui::SetNextWindowPos(ImVec2(20, display_h-fontSize*(1.1*numrows+1)));
      ImGuiWindowFlags window_flags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove;
      ImGui::Begin("Statistics", &show_stats_window, window_flags);
      ImGui::Text("Step %13ld", sim.get_nstep());
      ImGui::Text("Time %13.4f", sim.get_time());
      if (sim.get_npanels() > 0) ImGui::Text("Panels %11ld", sim.get_npanels());
      ImGui::Text("Particles %8ld", sim.get_nparts());
      if (sim.get_nfldpts() > 0) ImGui::Text("Field Pts %8ld", sim.get_nfldpts());
      ImGui::End();
    }

    // Show the terminal output of the program
    if (show_terminal_window) {
      ImGui::SetNextWindowSize(ImVec2(200,100), ImGuiSetCond_FirstUseEver);
      ImGui::Begin("Terminal", &show_terminal_window);
      ImGui::Text("Hello");
      ImGui::End();
    }

    // Show the ImGui test window. Most of the sample code is in ImGui::ShowTestWindow()
    if (show_test_window) {
      ImGui::SetNextWindowPos(ImVec2(650, 20), ImGuiSetCond_FirstUseEver);
      ImGui::ShowTestWindow();
    }

    // draw the simulation: panels and particles
    compute_ortho_proj_mat(window, rparams.vcx, rparams.vcy, &rparams.vsize, gl_projection);
    sim.drawGL(gl_projection, rparams);

    // if simulation has not been initted, draw the features instead!
    //for (auto const& bf : bfeatures) { bf.drawGL(gl_projection, rparams); }

    // here is where we write the buffer to a file
    if ((is_ready and record_all_frames and sim_is_running) or draw_this_frame) {
      static int frameno = 0;
      std::stringstream pngfn;
      pngfn << "img_" << std::setfill('0') << std::setw(5) << frameno << ".png";
      png_out_file = pngfn.str();
      (void) saveFramePNG(png_out_file);
      std::cout << "Wrote screenshot to " << png_out_file << std::endl;
      frameno++;
      draw_this_frame = false;
    }

    // if we're just drawing this one frame, then announce that we wrote it
    if (not png_out_file.empty()) {
      static int32_t pngframect = 0;
      ++pngframect;

      // draw the notification
      ImGui::SetNextWindowSize(ImVec2(10+fontSize*12, 10+fontSize*2));
      ImGui::SetNextWindowPos(ImVec2(0,0), 0);
      ImGuiWindowFlags window_flags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoResize;
      ImGui::Begin("Png written", NULL, window_flags);
      ImGui::Text("Wrote %s", png_out_file.c_str());
      ImGui::End();

      // make sure this isn't up for too long
      if (pngframect == 90) {
        pngframect = 0;
        png_out_file.clear();
      }
    }

    // draw the GUI
    ImGui::Render();
    // int display_w;
    // int display_h;
    // glfwMakeContextCurrent(window);
    // glfwGetFrameBufferSize(window, &display_w, &display_h);
    // glViewport(0, 0, display_w, display_h);
    // glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
    // glClear(GL_COLOR_BUFFER_BIT);
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
    // all done! swap buffers to the user can see
    glfwMakeContextCurrent(window);
    glfwSwapBuffers(window);
  }

  // Cleanup
  std::cout << "Starting shutdown procedure" << std::endl;
  sim.reset();
  std::cout << "Quitting" << std::endl;
  ImGui_ImplOpenGL3_Shutdown();
  ImGui_ImplGlfw_Shutdown();
  ImGui::DestroyContext();
  glfwDestroyWindow(window);
  glfwTerminate();

  return 0;
}

