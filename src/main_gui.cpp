/*
 * main_gui.cpp - Driver code for Omega2D + ImGui + Vc vortex particle method
 *                and boundary element method solver, GUI version
 *
 * (c)2017-20 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
 *            Blake B Hillier <blakehillier@mac.com>
 */

#include "FlowFeature.h"
#include "BoundaryFeature.h"
#include "MeasureFeature.h"
#include "Simulation.h"
#include "JsonHelper.h"
#include "Body.h"
#include "RenderParams.h"
#include "FeatureDraw.h"
#include "json/json.hpp"
#include "main_gui_functions.cpp"

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


int main(int argc, char const *argv[]) {
  std::cout << std::endl << "Omega2D GUI" << std::endl;

  // Set up vortex particle simulation
  Simulation sim;
  std::vector< std::unique_ptr<FlowFeature> > ffeatures;
  std::vector< std::unique_ptr<BoundaryFeature> > bfeatures;
  std::vector< std::unique_ptr<MeasureFeature> > mfeatures;
  FeatureDraw bdraw;
  size_t nframes = 0;
  static bool sim_is_running = false;
  static bool begin_single_step = false;

  // placeholder for command-line input file
  std::string command_line_input;

  // Set up primary OpenGL window
  glfwSetErrorCallback(error_callback);
  bool init = glfwInit();
  if (!init) {
    std::cout << "glfwInit failed" << std::endl;
    exit(-1);
  }

#if __APPLE__
  const char* glsl_version = "#version 150";
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#elif _WIN32
  const char* glsl_version = "#version 330 core";
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#else
  const char* glsl_version = "#version 150";
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#endif
  GLFWwindow* window = glfwCreateWindow(1280, 720, "Omega2D GUI", nullptr, nullptr);
  if (!window) { 
  std::cout << "glfwCreateWindow created " << window << std::endl;
  exit(-1);
  }
  glfwMakeContextCurrent(window);
  glfwSwapInterval(1); // Enable vsync

  //gl3wInit();
  init = gladLoadGL();
  if (!init) {
    std::cout << "gladLoadGL failed" << std::endl;
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
  bool save_all_bdry = false;		// save Boundary Features every step
  bool save_all_flow = false;		// save Flow Features every step
  bool save_all_meas = false;		// save Measure Features every step
  bool save_all_vtus = false;		// save all collections the coming step
  bool export_vtk_this_frame = false;	// write set of vtu files with the current data
  std::vector<std::string> vtk_out_files; // list of just-output files
  bool save_all_imgs = false;		// save screenshot every step
  bool export_png_when_ready = false;	// write frame to png as soon as its done
  bool write_png_immediately = false;	// write frame to png right now
  std::string png_out_file;		// the name of the recently-written png
  bool show_stats_window = true;
  bool show_welcome_window = true;
  bool show_terminal_window = false;
  bool show_demo_window = false;
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

  // Load file names and paths of pre-stored sims
  std::vector<nlohmann::json> sims;
  std::vector<std::string> descriptions = {"Select a simulation"};
  LoadJsonSims(sims, descriptions, EXAMPLES_DIR);

  // Main loop
  std::cout << "Starting main loop" << std::endl;
  while (!glfwWindowShouldClose(window)) {
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
        // and reset
        sim.reset();
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

      // split on which to write
      if (save_all_vtus) {
        // default is to save all collections to vtu files
        vtk_out_files = sim.write_vtk();
        // and don't do this next time
        save_all_vtus = false;

      } else if (save_all_bdry or save_all_flow or save_all_meas) {
        // only write select vtu files and don't echo
        (void) sim.write_vtk(-1, save_all_bdry, save_all_flow, save_all_meas);
      }

      // tell this routine next time around not to print
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

      if (save_all_bdry or save_all_flow or save_all_meas) export_vtk_this_frame = true;
      if (save_all_imgs) export_png_when_ready = true;

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
      ImGui::SetNextWindowSize(ImVec2(400,200), ImGuiCond_FirstUseEver);
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

    ImGui::SetNextWindowSize(ImVec2(140+fontSize*24,100+fontSize*12), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowPos(ImVec2(20,20), ImGuiCond_FirstUseEver);
    ImGui::Begin("Omega2D");
    ImGui::Spacing();

    // Select pre-populated simulations
    {
      int currentItemIndex = 0;
      const char* currentItem = descriptions[currentItemIndex].c_str();
      static ImGuiComboFlags flags = 0;
      if (ImGui::BeginCombo("", currentItem, flags)) // The second parameter is the label previewed before opening the combo.
      {
        for (size_t n = 0; n < descriptions.size(); n++)
        {
          bool is_selected = (currentItem == descriptions[n].c_str());
          if (ImGui::Selectable(descriptions[n].c_str(), is_selected)) {
            currentItem = descriptions[n].c_str();
            currentItemIndex = n;
          }
          if (is_selected) {
              ImGui::SetItemDefaultFocus();   // Set the initial focus when opening the combo (scrolling + for keyboard navigation support in the upcoming navigation branch)
          }
        }
        ImGui::EndCombo();
      }

      if (currentItemIndex) {
        sim.reset();
        bfeatures.clear();
        ffeatures.clear();
        mfeatures.clear();
        parse_json(sim, ffeatures, bfeatures, mfeatures, rparams, sims[currentItemIndex-1]);
        // clear and remake the draw geometry
        bdraw.clear_elements();
        for (auto const& bf : bfeatures) {
          bdraw.add_elements( bf->get_draw_packet(), bf->is_enabled() );
        }
        // finish setting up and run
        is_viscous = sim.get_diffuse();
        currentItemIndex = 0;
        sim_is_running = true;
      }
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
          nlohmann::json j = read_json(infile);
          parse_json(sim, ffeatures, bfeatures, mfeatures, rparams, j);

          // clear and remake the draw geometry
          bdraw.clear_elements();
          for (auto const& bf : bfeatures) {
            bdraw.add_elements( bf->get_draw_packet(), bf->is_enabled() );
          }

          // finish setting up and run
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
      nlohmann::json j = read_json(command_line_input);
      parse_json(sim, ffeatures, bfeatures, mfeatures, rparams, j);

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
    } // end Simulation Globals

    ImGui::Spacing();
    //if (ImGui::CollapsingHeader("Flow structures", ImGuiTreeNodeFlags_DefaultOpen)) 
    if (ImGui::CollapsingHeader("Startup structures")) {

      if (ffeatures.size() + bfeatures.size() == 0) {
        ImGui::Text("Add flow or boundry features (like vortex blobs and solid objects) here, then click RUN.");
      }


      ImGui::Spacing();

      // button and modal window for adding new boundary objects
      if (ImGui::Button("Add boundary")) ImGui::OpenPopup("New boundary structure");
      ImGui::SetNextWindowSize(ImVec2(400,275), ImGuiCond_FirstUseEver);
      if (ImGui::BeginPopupModal("New boundary structure"))
      {
        // define movement first
        static int mitem = 0;
        static char strx[512] = "0.0*t";
        static char stry[512] = "0.0*t";
        static char strrad[512] = "0.0*t";
        int changed = obj_movement_gui(mitem, strx, stry, strrad);

        // define geometry second
        static int item = 0;
        static int numItems = 7;
        const char* items[] = { "circle", "square", "oval", "rectangle", "segment", "polygon", "NACA 4-digit" };
        ImGui::Spacing();
        ImGui::Combo("geometry type", &item, items, numItems);

        // static bp prevents a bunch of pointers from being created during the same boundary creation
        // The switch prevents constant assignment (mainly to prevent the terminal from being flooded from messages)
        static std::shared_ptr<Body> bp = nullptr;
        if (changed) {
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
               break;
          }
        }

        if (BoundaryFeature::draw_creation_gui(item, mitem, sim.get_ips(), bp, bfeatures)) {
          if (mitem == 2) { sim.add_body(bp); }
          bfeatures.back()->generate_draw_geom();
          bdraw.add_elements( bfeatures.back()->get_draw_packet(), bfeatures.back()->is_enabled() );
        }
        ImGui::EndPopup();
      } // end new boundary structure

      // button and modal window for adding new flow structures
      ImGui::SameLine();
      if (ImGui::Button("Add vortex")) ImGui::OpenPopup("New flow structure");
      ImGui::SetNextWindowSize(ImVec2(400,200), ImGuiCond_FirstUseEver);
      if (ImGui::BeginPopupModal("New flow structure"))
      {
        static int item = 1;
        const char* items[] = { "single particle", "round vortex blob", "Gaussian vortex blob", "asymmetric vortex blob", "block of vorticity", "random particles", "particle emitter" };
        ImGui::Combo("type", &item, items, 7);

        // show different inputs based on what is selected
        switch(item) {
          case 0: {
            // creates a single particle
            SingleParticle::draw_creation_gui(ffeatures);
          } break;
         case 1: {
              // creates a blob of vorticies
              VortexBlob::draw_creation_gui(ffeatures, sim.get_ips());
          } break;
          case 2: {
            // a gaussian blob of multiple vorticies
            GaussianBlob::draw_creation_gui(ffeatures, sim.get_ips());
          } break;
          case 3: {
            // an asymmetric blob of multiple vorticies
            AsymmetricBlob::draw_creation_gui(ffeatures, sim.get_ips());
          } break;
          case 4: {
            // particles in a rectangle
            UniformBlock::draw_creation_gui(ffeatures, sim.get_ips());
          } break;
          case 5: {
            // random particles in a rectangle
            BlockOfRandom::draw_creation_gui(ffeatures);
          } break;
          case 6: {
            // create a particle emitter
            ParticleEmitter::draw_creation_gui(ffeatures);
          } break;
        }

        if (ImGui::Button("Cancel", ImVec2(120,0))) { ImGui::CloseCurrentPopup(); }
        ImGui::EndPopup();
      } // end popup new flow structures


      // button and modal window for adding new measurement objects
      ImGui::SameLine();
      if (ImGui::Button("Add measurement")) ImGui::OpenPopup("New measurement structure");
      ImGui::SetNextWindowSize(ImVec2(400,200), ImGuiCond_FirstUseEver);
      if (ImGui::BeginPopupModal("New measurement structure"))
      {
        static int item = 0;
        const char* items[] = { "single point/tracer", "streakline", "circle of tracers", "line of tracers", "measurement line", "measurement grid" };
        ImGui::Combo("type", &item, items, 6);

        // show different inputs based on what is selected
        switch(item) {
          case 0: {
            // a single measurement point
            SinglePoint::draw_creation_gui(mfeatures);
          } break;
          case 1: {
            // a tracer emitter
            TracerEmitter::draw_creation_gui(mfeatures);
          } break;
          case 2: {
            // a tracer circle
            TracerBlob::draw_creation_gui(mfeatures, rparams.tracer_scale, sim.get_ips());
          } break;
          case 3: {
            // a tracer line
            TracerLine::draw_creation_gui(mfeatures, rparams.tracer_scale, sim.get_ips());
          } break;
          case 4: {
            // a static, measurement line
            MeasurementLine::draw_creation_gui(mfeatures, rparams.tracer_scale, sim.get_ips());
          } break;
          case 5: {
            // a static grid of measurement points
            GridPoints::draw_creation_gui(mfeatures, sim.get_ips());
          } break;
        }

        if (ImGui::Button("Cancel", ImVec2(120,0))) { ImGui::CloseCurrentPopup(); }
        ImGui::EndPopup();
      } // end measurement structures 

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
        const bool ischeck = ImGui::Checkbox("", bfeatures[i]->addr_enabled());
        ImGui::PopID();
        if (bfeatures[i]->is_enabled()) {
          ImGui::SameLine();
          ImGui::Text("%s", bfeatures[i]->to_string().c_str());
        } else {
          ImGui::SameLine();
          ImGui::TextColored(ImVec4(0.5f,0.5f,0.5f,1.0f), "%s", bfeatures[i]->to_string().c_str());
        }

        // if the checkbox flipped positions this frame, ischeck is 1
        if (ischeck) bdraw.reset_enabled(i,bfeatures[i]->is_enabled());

        // add a "remove" button at the end of the line (so it's not easy to accidentally hit)
        ImGui::SameLine();
        ImGui::PushID(++buttonIDs);
        if (ImGui::SmallButton("remove")) del_this_bdry = i;
        ImGui::PopID();
      }
      if (del_this_bdry > -1) {
        std::cout << "Asked to delete boundary feature " << del_this_bdry << std::endl;
        bfeatures.erase(bfeatures.begin()+del_this_bdry);

        // clear out and re-make all boundary draw geometry
        bdraw.clear_elements();
        for (auto const& bf : bfeatures) {
          bdraw.add_elements( bf->get_draw_packet(), bf->is_enabled() );
        }
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
    if (ImGui::CollapsingHeader("Rendering controls")) { draw_render_gui(rparams); }
    
    // Solver parameters, under its own header
    ImGui::Spacing();
    if (ImGui::CollapsingHeader("Solver parameters (advanced)")) { sim.draw_advanced(); }

    // Output buttons, under a header
    ImGui::Spacing();
    if (ImGui::CollapsingHeader("Save output")) {

      // save setup
      ImGui::Spacing();
      ImGui::Text("Save simulation setup:");
      ImGui::SameLine();
      if (ImGui::Button("to JSON", ImVec2(10+4*fontSize,0))) show_file_output_window = true;

      // save current data
      ImGui::Separator();
      ImGui::Spacing();
      ImGui::Text("Save current step:");

      ImGui::SameLine();
      if (ImGui::Button("All to VTU", ImVec2(10+7*fontSize,0))) {
        save_all_vtus = true;
        export_vtk_this_frame = true;
      }
      ImGui::SameLine();
      if (ImGui::Button("Screenshot to PNG", ImVec2(10+10*fontSize,0))) write_png_immediately = true;

      // save data regularly
      ImGui::Separator();
      ImGui::Spacing();
      ImGui::Text("Save every step:");

      ImGui::Indent();
      ImGui::Checkbox("Boundary features (to VTU)", &save_all_bdry);
      ImGui::Checkbox("Flow features (to VTU)", &save_all_flow);
      ImGui::Checkbox("Measure features (to VTU)", &save_all_meas);
      ImGui::Checkbox("Window screenshot (to PNG)", &save_all_imgs);
      ImGui::Unindent();
    }

    if (show_file_output_window) {
      bool try_it = false;
      static std::string outfile = "file_name.json";

      if (fileIOWindow( try_it, outfile, recent_json_files, "Save", {"*.json"}, false, ImVec2(200+26*fontSize,300))) {
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
          write_png_immediately = true;

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

      if (ImGui::Button("ImGui Samples")) show_demo_window ^= 1;
      // use ASCII table for number: http://www.asciitable.com/
      // but use CAPITAL letter for a letter, jesus, really?!?
      //if (ImGui::IsKeyPressed(84) and not show_file_output_window) show_demo_window ^= 1;

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
    if (show_welcome_window) { show_welcome_window = draw_welcome_window(io.DisplaySize.x, io.DisplaySize.y); }

    // Show the simulation stats in the corner
    //if (nframes < 10){ std::cout << "show_stats_window: " << show_stats_window << std::endl; }
    if (show_stats_window) { draw_stats_window(sim.get_npanels(), sim.get_nfldpts(), sim.get_nstep(), 
                                               sim.get_time(), sim.get_nparts(), &show_stats_window,
                                               fontSize, display_h); }

    // Show the terminal output of the program
    if (show_terminal_window) {
      ImGui::SetNextWindowSize(ImVec2(200,100), ImGuiCond_FirstUseEver);
      ImGui::Begin("Terminal", &show_terminal_window);
      ImGui::Text("Hello");
      ImGui::End();
    }

    // Show the ImGui test window. Most of the sample code is in ImGui::ShowTestWindow()
    if (show_demo_window) {
      ImGui::SetNextWindowPos(ImVec2(650, 20), ImGuiCond_FirstUseEver);
      ImGui::ShowDemoWindow();
    }

    // draw the simulation: panels and particles
    compute_ortho_proj_mat(window, rparams.vcx, rparams.vcy, &rparams.vsize, gl_projection);
    sim.drawGL(gl_projection, rparams);

    // if simulation has not been initted, draw the features instead!
    if (not sim.is_initialized()) {
      // append draw geometries to FeatureDraw object
      for (auto const& bf : bfeatures) {
        if (bf->is_enabled()) {
          // what should we do differently?
        }
      }

      // and draw
      bdraw.drawGL(gl_projection, rparams);
    }

    // here is where we write the buffer to a file
    if ((is_ready and export_png_when_ready) or write_png_immediately) {
      static int frameno = 0;
      std::stringstream pngfn;
      pngfn << "img_" << std::setfill('0') << std::setw(5) << frameno << ".png";
      png_out_file = pngfn.str();
      (void) saveFramePNG(png_out_file);
      std::cout << "Wrote screenshot to " << png_out_file << std::endl;
      frameno++;
      // no need to tell the user every frame
      if (export_png_when_ready) png_out_file.clear();
      write_png_immediately = false;
      export_png_when_ready = false;
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

