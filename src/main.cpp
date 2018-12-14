/*
 * main.cpp - Driver code for Omega2D + ImGui + Vc vortex particle method and boundary element method solver
 *
 * (c)2017-8 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#include "FlowFeature.h"
#include "BoundaryFeature.h"
#include "Simulation.h"

#ifdef _WIN32
  // for glad
  #define APIENTRY __stdcall
  // for C++11 stuff
  #include <ciso646>
#endif
#include "glad.h"

#include "imgui/imgui.h"
#include "imgui/imgui_impl_glfw_gl3.h"

//#include <GL/gl3w.h>    // This example is using gl3w to access OpenGL
// functions (because it is small). You may use glew/glad/glLoadGen/etc.
// whatever already works for you.
#include <GLFW/glfw3.h>

#include <cmath>
#include <cstdio>
#include <iostream>
#include <memory>
#include <vector>

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
    (*_size) *= pow(1.05f, io.MouseWheel);

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
void compute_projection_matrix(GLFWwindow*         _thiswin,
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

int main(int, char**)
{
  // Set up vortex particle simulation
  Simulation sim;
  std::vector< std::unique_ptr<FlowFeature> > ffeatures;
  std::vector< std::unique_ptr<BoundaryFeature> > bfeatures;
  static bool sim_is_running = false;
  static bool begin_single_step = false;

  // Set up primary OpenGL window
  glfwSetErrorCallback(error_callback);
  if (!glfwInit())
    return 1;
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#if __APPLE__
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif
  GLFWwindow* window = glfwCreateWindow(1280, 720, "Omega2D with ImGui", nullptr, nullptr);
  glfwMakeContextCurrent(window);
  glfwSwapInterval(1); // Enable vsync

  //gl3wInit();

  if (!gladLoadGL()) {
    printf("Something went wrong!\n");
    exit(-1);
  }

  // Setup ImGui binding
  ImGui_ImplGlfwGL3_Init(window, true);

  //glfwSetKeyCallback(keyboard_callback);

  // Load Fonts
  // (there is a default font, this is only if you want to change it. see extra_fonts/README.txt for more details)

  // Get and set some IO functions
  ImGuiIO& io = ImGui::GetIO();
  io.IniFilename = ".omega2d.ini";

  // a string to hold any error messages
  std::string sim_err_msg;

  // projection matrix for the particles
  float vcx = -0.5f;
  float vcy = 0.0f;
  float vsize = 2.0f;
  std::vector<float> gl_projection;
  compute_projection_matrix(window, vcx, vcy, &vsize, gl_projection);

  // GUI and drawing parameters
  bool show_stats_window = false;
  bool show_terminal_window = false;
  bool show_test_window = false;
  ImVec4 clear_color = ImColor(15, 15, 15);
  ImVec4 pos_circ_color = ImColor(207, 47, 47);
  ImVec4 neg_circ_color = ImColor(63, 63, 255);
  //static bool show_origin = true;
  static bool is_viscous = false;

  // Main loop
  while (!glfwWindowShouldClose(window))
  {
    glfwPollEvents();
    ImGui_ImplGlfwGL3_NewFrame();

    //
    // Update simulation
    //

    // get results of latest step, if it just completed
    bool is_ready = sim.test_for_new_results();

    // see if we should start a new step
    //std::cout << "is_ready " << is_ready << "  sim_is_running " << sim_is_running << "  begin_single_step " << begin_single_step << std::endl;

    if (is_ready and (sim_is_running || begin_single_step)) {

      // if particles are not yet created, make them
      if (not sim.is_initialized()) {
        std::cout << "Initializing simulation" << std::endl;

        // initialize particle distributions
        for (auto const& ff: ffeatures) {
          sim.add_particles( ff->init_particles(sim.get_ips()) );
        }

        // initialize solid objects
        for (auto const& bf : bfeatures) {
          sim.add_boundary( bf->get_type(), bf->get_params() );
        }

        // initialize panels
        sim.init_bcs();
        sim.set_initialized();
      }

      // check flow for blow-up or errors
      sim_err_msg = sim.check_simulation(ffeatures.size(), bfeatures.size());

      if (sim_err_msg.empty()) {
        // the last simulation step was fine, OK to continue

        // generate new particles from emitters
        for (auto const& ff : ffeatures) {
          sim.add_particles( ff->step_particles(sim.get_ips()) );
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
      mouse_callback(window, io, &vcx, &vcy, &vsize);
    }

    // check for keypresses to toggle state
    //if (not io.WantCaptureKeyboard) {
      //keyboard_callback(
    //}

    //
    // The main Omega2D window
    //
    {

      ImGui::Begin("Omega2D");
      ImGui::TextWrapped("Welcome to Omega2D. First set your simulation globals, then add one or more flow structures, then press RUN. Have fun!");

      ImGui::Separator();
      ImGui::Text("Simulation globals");
      ImGui::SliderFloat("Time step", sim.addr_dt(), 0.001f, 0.1f, "%.4f", 2.0f);
      ImGui::Checkbox("Fluid is viscous (diffuses)", &is_viscous);
      if (is_viscous) {
        // show the toggle for AMR
        //static bool use_amr = false;
        //ImGui::Checkbox("Allow adaptive resolution", &use_amr);
        //sim.set_amr(use_amr);
        sim.set_diffuse(true);
        // and let user choose Reynolds number
        ImGui::SliderFloat("Reynolds number", sim.addr_re(), 10.0f, 10000.0f, "%.1f", 2.0f);
        ImGui::Text("Particle spacing %g", sim.get_ips());
      } else {
        static float my_ips = 0.03141;
        ImGui::SliderFloat("Particle spacing", &my_ips, 0.001f, 0.1f, "%.3f", 2.0f);
        // change underlying Re when this changes
        sim.set_re_for_ips(my_ips);
      }
      ImGui::InputFloat2("Freestream speed", sim.addr_fs());

      ImGui::Separator();
      ImGui::Text("Flow structures");

      int buttonIDs = 10;

      // list existing flow features here
      int del_this_item = -1;
      for (int i=0; i<(int)ffeatures.size(); ++i) {
        // add a "remove" button here somehow
        ImGui::PushID(++buttonIDs);
        if (ImGui::SmallButton("remove")) del_this_item = i;
        ImGui::PopID();

        //ImGui::SameLine();
        //ImGui::PushID(++buttonIDs);
        //if (ImGui::SmallButton("edit", ImVec2(60,0))) del_this_item = 0;
        //ImGui::PopID();

        ImGui::SameLine();
        ImGui::Text("%s", ffeatures[i]->to_string().c_str());
      }
      if (del_this_item > -1) {
        std::cout << "Asked to delete flow feature " << del_this_item << std::endl;
        ffeatures.erase(ffeatures.begin()+del_this_item);
      }

      // list existing boundary features here
      int del_this_bdry = -1;
      for (int i=0; i<(int)bfeatures.size(); ++i) {
        // add a "remove" button here somehow
        ImGui::PushID(++buttonIDs);
        if (ImGui::SmallButton("remove")) del_this_bdry = i;
        ImGui::PopID();

        ImGui::SameLine();
        ImGui::Text("%s", bfeatures[i]->to_string().c_str());
      }
      if (del_this_bdry > -1) {
        std::cout << "Asked to delete boundary feature " << del_this_bdry << std::endl;
        bfeatures.erase(bfeatures.begin()+del_this_bdry);
      }

      // button and modal window for adding new ones
      if (ImGui::Button("Add new flow structure")) ImGui::OpenPopup("New flow structure");
      ImGui::SetNextWindowSize(ImVec2(400,200), ImGuiSetCond_FirstUseEver);
      if (ImGui::BeginPopupModal("New flow structure"))
      {
        static int item = 1;
        const char* items[] = { "single particle", "vortex blob", "random particles", "particle emitter" };
        ImGui::Combo("type", &item, items, 4);

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
          case 0:
            // a single vortex particle
            ImGui::SliderFloat("strength", &str, -1.0f, 1.0f, "%.4f");
            ImGui::TextWrapped("This feature will add 1 particle");
            if (ImGui::Button("Add single particle")) {
              // this is C++14
              ffeatures.emplace_back(std::make_unique<SingleParticle>(xc[0], xc[1], str));
              // this is C++11
              //ffeatures.emplace_back(std::unique_ptr<SingleParticle>(new SingleParticle(xc[0], xc[1], str)));
              std::cout << "Added " << (*ffeatures.back()) << std::endl;
              ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            // it would be nice to be able to put this all in
            //SingleParticle::draw_creation_gui();
            break;

          case 1:
            // a blob of multiple vorticies
            ImGui::SliderFloat("strength", &str, -5.0f, 5.0f, "%.4f");
            ImGui::SliderFloat("radius", &rad, sim.get_ips(), 1.0f, "%.4f");
            ImGui::SliderFloat("softness", &soft, sim.get_ips(), 1.0f, "%.4f");
            ImGui::TextWrapped("This feature will add about %d particles", (int)(0.785398175*pow((2*rad+soft)/sim.get_ips(), 2)));
            if (ImGui::Button("Add vortex blob")) {
              ffeatures.emplace_back(std::make_unique<VortexBlob>(xc[0], xc[1], str, rad, soft));
              //ffeatures.emplace_back(std::unique_ptr<VortexBlob>(new VortexBlob(xc[0], xc[1], str, rad, soft)));
              std::cout << "Added " << (*ffeatures.back()) << std::endl;
              ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            break;

          case 2:
            // random particles in a rectangle
            ImGui::SliderInt("number", &npart, 1, 10000);
            ImGui::SliderFloat2("box size", xs, 0.01f, 10.0f, "%.4f", 2.0f);
            ImGui::DragFloatRange2("strength range", &strlo, &strhi, 0.001f, -0.1f, 0.1f);
            ImGui::TextWrapped("This feature will add %d particles", npart);
            if (ImGui::Button("Add random vorticies")) {
              ffeatures.emplace_back(std::make_unique<BlockOfRandom>(xc[0], xc[1], xs[0], xs[1], strlo, strhi, npart));
              //ffeatures.emplace_back(std::unique_ptr<BlockOfRandom>(new BlockOfRandom(xc[0], xc[1], xs[0], xs[1], strlo, strhi, npart)));
              std::cout << "Added " << (*ffeatures.back()) << std::endl;
              ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            break;

          case 3:
            // create a particle emitter
            static float estr = 0.1f;
            ImGui::SliderFloat("strength", &estr, -0.1f, 0.1f, "%.4f");
            ImGui::TextWrapped("This feature will add 1 particle per time step");
            if (ImGui::Button("Add particle emitter")) {
              // this is C++11
              ffeatures.emplace_back(std::make_unique<ParticleEmitter>(xc[0], xc[1], estr));
              //ffeatures.emplace_back(std::unique_ptr<ParticleEmitter>(new ParticleEmitter(xc[0], xc[1], estr)));
              std::cout << "Added " << (*ffeatures.back()) << std::endl;
              ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            break;
        }

        if (ImGui::Button("Cancel", ImVec2(120,0))) { ImGui::CloseCurrentPopup(); }
        ImGui::EndPopup();
      }

      // button and modal window for adding new boundary objects
      ImGui::SameLine();
      if (ImGui::Button("Add new boundary structure")) ImGui::OpenPopup("New boundary structure");
      ImGui::SetNextWindowSize(ImVec2(400,200), ImGuiSetCond_FirstUseEver);
      if (ImGui::BeginPopupModal("New boundary structure"))
      {
        static int item = 0;
        //const char* items[] = { "solid circle", "solid square", "solid object from file", "draw outline in UI" };
        const char* items[] = { "solid circle", "solid square" };
        ImGui::Combo("type", &item, items, 2);

        static float xc[2] = {0.0f, 0.0f};

        // always ask for center
        ImGui::InputFloat2("center", xc);

        // show different inputs based on what is selected
        switch(item) {
          case 0:
            // create a circular boundary
            static float circdiam = 1.0;
            ImGui::SliderFloat("diameter", &circdiam, 0.01f, 10.0f, "%.4f", 2.0);
            ImGui::TextWrapped("This feature will add a solid circular body centered at the given coordinates");
            if (ImGui::Button("Add circular body")) {
              bfeatures.emplace_back(std::make_unique<SolidCircle>(xc[0], xc[1], circdiam));
              //bfeatures.emplace_back(std::unique_ptr<SolidCircle>(new SolidCircle(xc[0], xc[1], circdiam)));
              std::cout << "Added " << (*bfeatures.back()) << std::endl;
              ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            break;
          case 1:
            // create a square/rectangle boundary
            static float sqside = 1.0;
            ImGui::SliderFloat("side length", &sqside, 0.1f, 10.0f, "%.4f");
            ImGui::TextWrapped("This feature will add a solid square body centered at the given coordinates");
            if (ImGui::Button("Add square body")) {
              ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            break;
/*
          case 2:
            // load a solid object as an outline
            ImGui::TextWrapped("This feature is currently unsupported. Sorry!");
            ImGui::SameLine();
            break;
*/
        }

        if (ImGui::Button("Cancel", ImVec2(120,0))) { ImGui::CloseCurrentPopup(); }
        ImGui::EndPopup();
      }

      ImGui::Separator();
      ImGui::Text("Rendering parameters");
      ImGui::ColorEdit3("positive circulation", (float*)&pos_circ_color);
      ImGui::ColorEdit3("negative circulation", (float*)&neg_circ_color);
      ImGui::ColorEdit3("background color", (float*)&clear_color);
      //ImGui::Checkbox("show origin", &show_origin);

      if (ImGui::Button("Recenter")) {
        // put everything back to center
        vcx = -0.5f;
        vcy = 0.0f;
        vsize = 2.0f;
      }


      ImGui::Separator();
      if (sim_is_running) {
        ImGui::Text("Simulation is running...time = %g", sim.get_time());
        if (ImGui::Button("PAUSE", ImVec2(200,0))) sim_is_running = false;
      } else {
        ImGui::Text("Simulation is not running, time = %g", sim.get_time());
        if (ImGui::Button("RUN", ImVec2(200,0))) sim_is_running = true;
        ImGui::SameLine();
        if (ImGui::Button("Step", ImVec2(120,0))) begin_single_step = true;
      }
      ImGui::SameLine();
      if (ImGui::Button("Reset", ImVec2(120,0))) {
        std::cout << std::endl << "Reset requested" << std::endl;
        // remove all particles and reset timer
        sim.reset();
        std::cout << "Reset complete" << std::endl;
      }

      ImGui::Separator();
      /*
      ImGui::Text("Open additional windows");
      if (ImGui::Button("Plot statistics")) show_stats_window ^= 1;
      ImGui::SameLine();
      if (ImGui::Button("Show terminal output")) show_terminal_window ^= 1;
      ImGui::SameLine();
      if (ImGui::Button("ImGui Samples")) show_test_window ^= 1;
      //if (ImGui::Button("Another Window")) show_another_window ^= 1;
      */
      ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);

      ImGui::Text("Number of panels: %ld  Number of particles: %ld", sim.get_npanels(), sim.get_nparts());

      ImGui::End();
    }

    // Show the simulation stats as 2D plots
    if (show_stats_window)
    {
      ImGui::SetNextWindowSize(ImVec2(200,100), ImGuiSetCond_FirstUseEver);
      ImGui::Begin("Statistics", &show_stats_window);
      ImGui::Text("Hello");
      ImGui::End();
    }

    // Show the terminal output of the program
    if (show_terminal_window)
    {
      ImGui::SetNextWindowSize(ImVec2(200,100), ImGuiSetCond_FirstUseEver);
      ImGui::Begin("Terminal", &show_terminal_window);
      ImGui::Text("Hello");
      ImGui::End();
    }

    // Show the ImGui test window. Most of the sample code is in ImGui::ShowTestWindow()
    if (show_test_window)
    {
      ImGui::SetNextWindowPos(ImVec2(650, 20), ImGuiSetCond_FirstUseEver);
      ImGui::ShowTestWindow(&show_test_window);
    }

    // Rendering
    int display_w, display_h;
    glfwGetFramebufferSize(window, &display_w, &display_h);
    glViewport(0, 0, display_w, display_h);
    glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
    glClear(GL_COLOR_BUFFER_BIT);

    // draw the simulation: panels and particles
    compute_projection_matrix(window, vcx, vcy, &vsize, gl_projection);
    sim.drawGL(gl_projection, &pos_circ_color.x, &neg_circ_color.x);

    // draw the GUI
    ImGui::Render();

    // all done! swap buffers to the user can see
    glfwSwapBuffers(window);
  }

  // Cleanup
  ImGui_ImplGlfwGL3_Shutdown();
  glfwTerminate();

  return 0;
}

