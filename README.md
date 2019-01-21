# Omega2D
Two-dimensional flow solver with GUI using vortex particle and boundary element methods

## Overview
[Computational Fluid Dynamics (CFD)](https://en.wikipedia.org/wiki/Computational_fluid_dynamics) encompasses a wide variety of methods to aid in the numerical simulation of fluid flows on digital computers. Most methods rely on the subdivision of the fluid domain into small, stationary cells, such as tetrahedra, and solve the [Navier-Stokes equations](https://en.wikipedia.org/wiki/Navier%E2%80%93Stokes_equations) on each Eulerian (not moving) cell. In contrast, vortex methods rely on a Lagrangian (moving with the flow) description of the only the [vorticity](https://en.wikipedia.org/wiki/Vorticity)-containing region of the fluid domain and any solid boundaries present. This eliminates many of the difficulties present in traditional CFD. In addition, the form of the equations used also removes the pressure term from the Navier-Stokes equations, which is a large source of instability and extra effort in traditional CFD. This is why many new flow solvers for unsteady momentum-dominated flows (non-microscopic in scale) are implemented using vortex methods.

Omega2D is a platform for testing methods and techniques for implementing a combined Lagrangian-Eulerian fluid flow solver for unsteady flows with complex boundaries. Lessons learned here will be applied to a full three-dimensional GUI and solver. The current version supports 2nd order convection and diffusion with uniform freestream and an arbitrary number of solid circular objects. The main focus right now is functionality instead of correctness.

This open-source code is aimed at users interested in understanding vortex methods as a tool for fluid simulation, or simply eager to try a fast fluid simulator without the gross approximations present in most other real-time tools.

## Build and run
This code uses some C++17 features, so should compile on GCC 7, Clang 4, and MSVC 19.10 compilers.

#### Prerequisites
Users will also need CMake and GLFW version 3 on their machines to build this, other requirements are included in this distribution. Get these on Fedora with

    sudo yum install cmake glfw3-devel eigen3-devel

or on Ubuntu with

    sudo apt-get install cmake glfw3-dev libeigen3-dev

or on OSX via [Homebrew](https://docs.brew.sh/Installation) with

    brew install cmake glfw eigen

#### Optional libraries
[Vc](https://github.com/VcDevel/Vc) is a vectorization library, and Omega2D uses it to greatly accelerate the velocity evaluations. This package can be built and installed external to Omega2D with

    git clone https://github.com/VcDevel/Vc.git
    cd Vc
    mkdir build
    cd build
    cmake -DCMAKE_INSTALL_PREFIX=/opt/Vc -DBUILD_TESTING=OFF ..
    make -j 4
    sudo make install
    cd ../..

#### Compile and run
Upon installation of the prerequisites, the following commands should build Omega2D.

    git clone git@github.com:Applied-Scientific-Research/Omega2D.git
    cd Omega2D
    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release -DUSE_OMP=ON -DUSE_VC=OFF ..
    make
    ./Omega2D.bin

If you were able to build and install Vc, then you should set `-DUSE_VC=ON` in the above `cmake` command.

On OSX, to get OpenMP parallization of the solver, you may need to use Clang, or install and use GCC with brew, and use a special `cmake` command:

    brew install gcc
    cmake -DCMAKE_C_COMPILER=/usr/local/bin/gcc-7 -DCMAKE_CXX_COMPILER=/usr/local/bin/g++-7 ..

## Set up and run a simulation

Upon running Omega2D, you will see the following GUI:

![screenshot](media/Screenshot_v4a.png?raw=true "Fresh GUI window")

In this window you can toggle viscous mode, set a freestream and time step, and tweak the simulation colors. But, most importantly, you can set the initial conditions, which consist of any number of *flow structures* and *boundary structures*. Pictured below is the action of adding a boundary structure:

![screenshot](media/Screenshot_v4b.png?raw=true "Adding a boundary structure")

Once flow and boundary structures are created, and the simulation parameters are set as desired, press the *RUN* button to start the simulation. At any time you can press *PAUSE* to pause the simulation or *Reset* to go back to the original conditions. At any time, you can left-click and drag on the flow field to move your point of view, or use the scroll wheel to zoom and unzoom.

Pictured below is a simulation of viscous flow over a circle at Reynolds number 200 after two steps. The blue and red fields represent negative and positive vorticity (rotation). This is created by flow sliding over the solid bondary. Because this flow solver uses vortex methods, we only require computational elements (vortex particles) where there is vorticity.

![screenshot](media/Screenshot_v4c.png?raw=true "2D flow over a circular cylinder at Re=200")

## To do
Tasks to consider or implement:

* Need to periodically check tracer particles to make sure they do not go inside of objects - like every frame we check 1/10th of all tracers and bump them out
* When background is white, nothing else shows up! Check blending mode. Need this to change so that we can make more attractive visuals. Like, presets for "technical (red/blue), b/w (white/grey), vibrant (??)
* If I add a measurement structure in the middle of a simulation, it doesn't init - should it?
* Add rendering parameters to the json file read and write, this includes window size and projection!
* Add file selector dialog, like [this one](https://github.com/gileoo/Imgui-IGS-Snippets) or [this one](https://github.com/Flix01/imgui/tree/imgui_with_addons/addons/imguifilesystem)
* Allow Points to be able to draw not only blobs but dots at the middles, too, using same arrays but different draw programs (turn either on or off?)
* Allow inert Points collections to never allocate space for radius
* Move some initialization back into ElementBase - like positions and such, keep radius in Points, then ElementBase can draw points?
* Move the GUI and pre-run parts of the various Body classes into their own class/file and out of main.cpp and Body.h
* When running, grey out the dt and Re fields - those are the only things you can't change
* Add a "ms/frame" and "FPS" for the simulation component also
* Add a batch mode also - once you get json reading and png writing working - do we really need this? yes.
* If we're adding png output, might as well add a "Record" button!
* Add "got it" button to first section (the welcome section) to make it go away (forever?)
* Add a "Run to time..." button to allow users to run a sim up to an exact time
* Add "Save setup", "Save flow", and "Save image" buttons
* Ideal initial interface: lots of stuff hidden, just a graphical menu with circles, squares, vortex patches, etc. Each has handles that you can drag to resize and reposition the element; all sizes/locations quantized to 0.1 or 0.05. "Expert" box lets you change Re, dt, etc.
* Add an animated GIF to the page to show how to set up a run?
* Have "status" line indicate when we're waiting for a step to finish after we've hit pause
* Instead of manipulating the projection matrix, have the mouse change the view matrix (assume model matrix is unity), see [here](https://solarianprogrammer.com/2013/05/22/opengl-101-matrices-projection-view-model/) for a nice write-up on the three OpenGL matrices
* Support faster VRM by caching values from first half to use in second half
* Make more of the sliders dynamic (like dt ~~and color~~) - be able to add new particles while a sim is running
* Right-click on the draw screen to add features - hard? can imgui handle it?
* Draw something when you add a feature (so we know it's doing something)
* Create an OpenGL compute shader routine for the particle-particle influence
* Add option to draw particles as thin white dots/lines "Draw elements"
* Add to Core.h a routine to precalculate and save the "trim tables" to be used when pushing particles away from boundaries
* Use these correct trim tables to push particles away from the boundary
* Re-orient the VRM insertion points to align to the nearest boundary - should smooth out the shedding
* Draw a freestream arrow in the LR corner
* Let the user grab the fs arrow to dynamically change the freestream
* Add field points in a grid over the visible domain, find vels, display as streaks
* Draw panels as polygons extending along the normal with a sharp edge on the body side and a gradient to zero on the flow side - then they can visually merge with the particles to make a visibly smooth and more-correct vorticity field
* Add Read Setup, Read State, Write Setup, Write State buttons, which call a JSON library, like [nlohmann/json](https://github.com/nlohmann/json), to perform those actions
* Reconsider templatizing on the scalar type. If you don't intend for floats or doubles in the same code, perhaps create a header with `using Scalar = float;` so you can flip back and forth easier. NBL
* Consider different method for including shader code NBL
* Consider reducing the number of virtual methods; the `variant` path NBL
* Recommend all private variables be prefixed with `m_` or suffixed with `_` to avoid name collisions with method parameters, etc, improve reading. NBL

Completed tasks:

* ~~Support multiple Lagrangian element collections under an Elements structure - new arch will allow this~~
* ~~Generalize bodies to allow squares and circles in the GUI~~
* ~~Add json reading and writing to allow batch modes - just like 3D code~~
* ~~Replace core architecture (using std::variant and Collection)~~
* ~~Get Vc into the influence calculations for the new arch~~
* ~~For new BEM: get particles to create correct RHS, then panels to create correct A matrix, then Points to generate correct particles, finally for Surfaces to remove Points beneath them~~
* ~~Get the Panels into the new arch - includes creation, BEM, vel-finding, drawing~~
* ~~Upgrade Merge to the new arch, then particle-only diffusion will be complete~~
* ~~Get new arch to perform VRM~~
* ~~Get new arch to create a line or a blob of inert particles, compute their motion, and draw them~~
* ~~Set tracer points to be a constant fraction of the minimum particle size~~
* ~~Echo correct number of particles and points in the GUI~~
* ~~Make alternative path for drawing points - do not pass strength or radius, just color and position~~
* ~~Make a second copy of all particles in Points and compute their motion, also - compare to Particles~~
* ~~Move radius from ElementBase to Points~~
* ~~Space bar pauses and resumes simulation (need to capture keyboard inputs)~~
* ~~Add warning when there are no flow structures or boundary structures~~
* ~~Add warning when there are bodies, but no freestream and no flow structures and no particles (because nothing will happen)~~
* ~~Make the sections of the GUI collapsable~~
* ~~Limit the upper level for dt to 0.1, lower level for Re at 10, limit circle radius between 0.01 and 10.0 but use logarithmic slider!~~
* ~~Find code that computes influence of thick-"cored" panels on a point - see Pate's dissertation - nope.~~
* ~~Move code to public repository~~

## Thanks
This project is funded by the [National Institutes of Health (NIH)](https://www.nih.gov/) under grant number 1 R01 EB022180-01A1 ("A Fast High-Order CFD for Turbulent Flow Simulation in Cardio-Devices").

Thanks to [Omar Cornut](http://www.miracleworld.net/) for his [dear imgui](https://github.com/ocornut/imgui) library, sol-prog's [OpenGL Tutorials](https://github.com/sol-prog/OpenGL-101), and Jim Susinno's [OpenGL-Boilerplate](https://github.com/jimbo00000/OpenGL-Boilerplate).

VRM code is functional thanks to jlblancoc for [Nanoflann](https://github.com/jlblancoc/nanoflann) (a header-only tree search library), and to all of the developers of [Eigen](http://eigen.tuxfamily.org/) (a C++ matrix/vector library). The BEM code also relies heavily on [Eigen](http://eigen.tuxfamily.org/). We also love [Vc](https://github.com/VcDevel/Vc), an excellent SIMD library by Matthias Kretz.

JSON reading and writing is thanks to [JSON for Modern C++](https://github.com/nlohmann/json) by [Niels Lohmann](http://nlohmann.me).
