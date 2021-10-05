# Omega2D
Two-dimensional flow solver with GUI using vortex particle and boundary element methods

*NOTE: This program currrently uses O(N^2) calculations for velocity, so runs more slowly than desired. A O(N log N) treecode is being prepared to alleviate this.*

![startupvideo](media/IntroCircle1.gif?raw=true "Session sample")

## Overview
[Computational Fluid Dynamics (CFD)](https://en.wikipedia.org/wiki/Computational_fluid_dynamics) encompasses a wide variety of methods to aid in the numerical simulation of fluid flows on digital computers. Most methods rely on the subdivision of the fluid domain into small, stationary cells, such as tetrahedra, and solve the [Navier-Stokes equations](https://en.wikipedia.org/wiki/Navier%E2%80%93Stokes_equations) on each Eulerian (not moving) cell. In contrast, vortex methods rely on a Lagrangian (moving with the flow) description of the only the [vorticity](https://en.wikipedia.org/wiki/Vorticity)-containing region of the fluid domain and any solid boundaries present. This eliminates many of the difficulties present in traditional CFD. In addition, the form of the equations used also removes the pressure term from the Navier-Stokes equations, which is a large source of instability and extra effort in traditional CFD. This is why many new flow solvers for unsteady momentum-dominated flows (non-microscopic in scale) are implemented using vortex methods.

Omega2D is a platform for testing methods and techniques for implementing a combined Lagrangian-Eulerian fluid flow solver for unsteady flows with complex boundaries. Lessons learned here will be applied to a full three-dimensional GUI and solver. The current version supports 2nd order convection and diffusion with uniform freestream, an arbitrary number of moving objects, and low to moderate Reynolds numbers. The main focus of effort right now is toward functionality instead of absolute correctness.

This open-source code is aimed at users interested in understanding vortex methods as a tool for fluid simulation, or simply eager to try a fast fluid simulator without the gross approximations present in most other real-time tools.


## Build the software
This code uses some C++17 features (like `std::variant` and `<filesystem>`, so requires GCC 8, Clang 5, and MSVC 19.10 (Visual Studio 15 2017) or newer compilers.

#### Prerequisites
Users will also need CMake, Eigen (version 3.3 or newer), and GLFW version 3 on their machines to build this, other requirements are included in this distribution. Get these on Fedora with

    sudo dnf install cmake glfw3-devel eigen3-devel

The hybrid version also requires the following packages.

    sudo dnf install libxml2-devel lapack-devel

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

The above commands should work verbatim on Linux and OSX. Don't ask about Windows - there's a calling convention issue preventing this from working.

#### Compile
Upon installation of the prerequisites, the following commands should build Omega2D.

    git clone --recurse-submodules git@github.com:Applied-Scientific-Research/Omega2D.git
    cd Omega2D
    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release -DUSE_HO_HYBRID=CXX -DUSE_OMP=ON -DUSE_VC=OFF ..
    make

If you were able to build and install Vc, then you should set `-DUSE_VC=ON` in the above `cmake` command.

To use the system Clang on Linux, you will want the following variables defined:

    cmake -DCMAKE_C_COMPILER=/usr/bin/clang -DCMAKE_CXX_COMPILER=/usr/bin/clang++ ..

On OSX, to get OpenMP parallelization of the solver, you may need to install GCC with brew (as above), and add a few more arguments to the `cmake` command:

    brew install gcc
    cmake -DCMAKE_C_COMPILER=/usr/local/bin/gcc-x -DCMAKE_CXX_COMPILER=/usr/local/bin/g++-x ..

where x is the latest version on your machine (you can check this by going to `/usr/local/bin`).

## Run a simulation in the GUI
If you were able to build the software, you should be able to run

    ./Omega2D.bin

Upon running Omega2D, you will see a GUI floating over a black field. Using the *Select a simulation...* pull-down menu, you can quickly load and run a preset simulation. Let's load "flow over circle".

![screenshot](media/Screenshot_v5a.png?raw=true "Load a preset simulation")

At any time you can press *PAUSE* to pause the simulation or *Reset* to go back to the original conditions. At any time, you can left-click and drag on the flow field to move your point of view, or use the scroll wheel to zoom and unzoom. Space bar also pauses and unpauses the simulation. Note that some simulations quickly become large enough to take several seconds between updates. Don't worry: when you pause, the current simulation step will finish.

There are several collapsible headers which you can open to modify this simulation, those include *Simulation globals* such as viscosity and flow speed, *Flow structures* such as solid bodies, vortex blobs, and tracers, and *Rendering parameters*. Some changes you make in these fields will affect the simulation immediately, but most will require you to *Reset*.

![screenshot](media/Screenshot_v5b.png?raw=true "Flow over a circular cylinder")

Pictured above is a simulation of viscous flow over a circular cylinder at Reynolds number 250 after 76 steps. The blue and red fields represent negative and positive vorticity (rotation). Vorticity is created when flow moves over a solid boundary, but must stick to the boundary surface. Because this flow solver uses vortex methods, we only require computational elements (vortex particles) where there is vorticity - nowhere else.

### Run a hybrid case
Omega2D now comes batched with HO-CXX, a submodule which contains a high-order Eulerian velocity-vorticity CFD solver. This method can be used to achieve better resolution and accuracy in regions with high vorticity gradients like walls.

To run a case which uses this Eulerian solver in some regions and the default Lagrangian vortex methods solver in all other regions, you will need a mesh created by a more recent version of [GMSH](https://gmsh.info/). Look in the `extern/gmsh-reader/data` directory for some scripts to assist you. This mesh must define one or more named `Physical Curve` objects called `open`, and optionally `wall`, `inlet`, and `outlet`. The fluid domain should be a `Physical Surface("fluid")`.

Load the mesh file into the GUI under `Startup structures`->`Add boundary`->`gmsh file`. Then under the `Solver parameters (advanced)` tab, inside `Hybrid/Grid settings` you can set runtime parameters such as the element order, integration order, and number of substeps.

Output from the gridded region will automatically be written to high-order xml-based VTK files when you click `All to VTU`. Open these files in [ParaView](https://www.paraview.org/) alongside the particle data.

### Run a batch job
If you already have an input file in JSON format, or you exported one from the GUI, you can run a batch (no GUI) simulation with

    ./Omega2Dbatch.bin input.json

Output will be written to the terminal and files to the working directory.

### Render a movie
The GUI has an option to `RECORD to png`. When you press this button, the simulation will progress as fast as it can, writing the flow field to a PNG image every time step. It is suggested that you set your view point first, then `Reset`, before recording.

Generate an X.264-encoded video from a series of png images with the following command. Make sure to use the actual resolution of the images. The reason for all the extra options is to ensure that the resulting video will play on Linux, Windows, and Mac - Quicktime is very picky about which video files it will play.

    mencoder "mf://img*png" -mf w=1280:h=720:type=png:fps=30 -o video.mp4 -sws 9 -of lavf -lavfopts format=mp4 -nosub -vf softskip,harddup -nosound -ovc x264 -x264encopts bitrate=4000:vbv_maxrate=6000:vbv_bufsize=2000:nointerlaced:force_cfr:frameref=3:mixed_refs:bframes=1:b_adapt=2:weightp=1:direct_pred=auto:aq_mode=1:me=umh:me_range=16:subq=6:mbtree:psy_rd=0.8,0.2:chroma_me:trellis=1:nocabac:deblock:partitions=p8x8,b8x8,i8x8,i4x4:nofast_pskip:nodct_decimate:threads=auto:ssim:psnr:keyint=300:keyint_min=30:level_idc=30:global_header

    ffmpeg -f image2 -pattern_type glob -i "img*png" -c:v libx264 -crf 20 -framerate 30 -profile:v baseline -level 3.0 -pix_fmt yuv420p -f mp4 video_264.mp4

    ffmpeg -f image2 -framerate 30 -pattern_type glob -i "img_*.png" -c:v libx265 -crf 25 -tune fastdecode -pix_fmt yuv420p -f mp4 video_265.mp4

## To do
Tasks to consider or implement:

* Move boundary conditions to their own class, with normal and tangential vels, inlets/outlets, constant/formula, etc.
* Add support for NACA wings to be created by reading geometry files in standard formats - this might mean enabling Kutta points (reactive Points)
* Add an openmp setting to the CPU execution environment, to show the effects of multithreading
* Bug: static field points should not be pushed out of bodies
* Compute vel grad at each point (vortex or field pt) and write it to the vtu file
* If we have velgrad, we can compute the elongation of any point - use this to determine when a field point has stretched too far, and replace it with two child field points, of half-brightness, in the correct places; this should be a toggleable behavior, as it would demand ever-increasing numbers of field points
* Have method to tell batch and gui to dump vtu files periodically
* Would be nice to have an option to have solid bodies randomly shed tracers (small chance per step per panel)
* Consider a zlib implementation for both vtu and png files: [zstr](https://github.com/mateidavid/zstr), [miniz](https://github.com/richgel999/miniz), or [gzip-hpp](https://github.com/mapbox/gzip-hpp)
* How awesome would it be to show the flow as a time-consistent LIC image? See UFLIC (Shen & Kao, IEEE ToVaCG 1999)
* When background is white, nothing else shows up! Check blending mode. Need this to change so that we can make more attractive visuals. Like, presets for "technical (red/blue), b/w (white/grey), vibrant (??)
* If I add a measurement structure in the middle of a simulation, it doesn't init - should it?
* Move some initialization back into ElementBase - like positions and such, keep radius in Points, then ElementBase can draw points?
* Add a "ms/frame" and "FPS" for the simulation component also
* Add "got it" button to first section (the welcome section) to make it go away (forever?)
* Ideal initial interface: lots of stuff hidden, just a graphical menu with circles, squares, vortex patches, etc. Each has handles that you can drag to resize and reposition the element; all sizes/locations quantized to 0.1 or 0.05. "Expert" box lets you change Re, dt, etc.
* Instead of manipulating the projection matrix, have the mouse change the view matrix (assume model matrix is unity), see [here](https://solarianprogrammer.com/2013/05/22/opengl-101-matrices-projection-view-model/) for a nice write-up on the three OpenGL matrices
* Support 2nd order time accuracy in VRM by caching values from first half to use in second half
* Make more of the sliders dynamic (like dt ~~and color~~) - be able to add new particles while a sim is running
* Right-click on the draw screen to add features - hard? can imgui handle it?
* Create an OpenGL compute shader routine for the particle-particle influence, like in Omega3D
* Add option to draw particles as thin white dots/lines "Draw elements"
* Draw a freestream arrow in the LR corner
* Let the user grab the fs arrow to dynamically change the freestream
* Allow a formula for the freestream instead of it being constant
* Add field points in a grid over the visible domain, find vels, display as streaks - make this one of a few rendering options
* Draw panels as polygons extending along the normal with a sharp edge on the body side and a gradient to zero on the flow side - then they can visually merge with the particles to make a visibly smooth and more-correct vorticity field
* Reconsider templatizing on the scalar type. If you don't intend for floats or doubles in the same code, perhaps create a header with `using Scalar = float;` so you can flip back and forth easier. NBL
* Consider different method for including shader code NBL
* Consider reducing the number of virtual methods; the `variant` path NBL
* Recommend all private variables be prefixed with `m_` or suffixed with `_` to avoid name collisions with method parameters, etc, improve reading. NBL


## Thanks
This project is funded by the [National Institutes of Health (NIH)](https://www.nih.gov/) under grant number 1 R01 EB022180-01A1 ("A Fast High-Order CFD for Turbulent Flow Simulation in Cardio-Devices").

Thanks to [Omar Cornut](http://www.miracleworld.net/) for his [dear imgui](https://github.com/ocornut/imgui) library, file browser dialogs from [Imgui-IGS-Snippets](https://github.com/gileoo/Imgui-IGS-Snippets), sol-prog's [OpenGL Tutorials](https://github.com/sol-prog/OpenGL-101), Jim Susinno's [OpenGL-Boilerplate](https://github.com/jimbo00000/OpenGL-Boilerplate), and the [miniz](https://github.com/richgel999/miniz/tree/master) compression library.

VRM code is functional thanks to jlblancoc for [Nanoflann](https://github.com/jlblancoc/nanoflann) (a header-only tree search library), and to all of the developers of [Eigen](http://eigen.tuxfamily.org/) (a C++ matrix/vector library). The BEM code also relies heavily on [Eigen](http://eigen.tuxfamily.org/). We also love [Vc](https://github.com/VcDevel/Vc), an excellent SIMD library by Matthias Kretz.

JSON reading and writing is thanks to [JSON for Modern C++](https://github.com/nlohmann/json) by [Niels Lohmann](http://nlohmann.me). XML output to VTK files is done using [tinyxml2](https://github.com/leethomason/tinyxml2) and [cppcodec](https://github.com/tplgy/cppcodec) for base64 encoding. And mathematical expression parsing came from [Lewis Van Winkle](https://codeplea.com/)'s [tinyexpr](https://github.com/codeplea/tinyexpr).

Many thanks to NBL for valuable discussions of architecture and C++ syntax and idioms.
