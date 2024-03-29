#
# Omega2D
# (c)2017-2023 Applied Scientific Research, Inc.
#
CMAKE_MINIMUM_REQUIRED( VERSION 3.4 )
PROJECT( Omega2D )

# User-serviceable variables here
SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo")
IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE "Release")
ENDIF ()
SET (CMAKE_BUILD_TYPE "Release" CACHE STRING "Select which configuration to build")
SET (BUILD_GUI TRUE CACHE BOOL "Build the GUI version")
SET (BUILD_BATCH FALSE CACHE BOOL "Build the batch (no GUI) version")
SET (CMAKE_INSTALL_PREFIX CACHE PATH "Installation location for binaries, sample inputs, and licenses")
SET (USE_OMP TRUE CACHE BOOL "Use OpenMP multithreading")
SET (USE_VC FALSE CACHE BOOL "Use Vc for vector arithmetic")
SET (USE_PLUGIN_AVRM FALSE CACHE BOOL "Enable experimental adaptive VRM plugin")
SET (USE_PLUGIN_SIMPLEX FALSE CACHE BOOL "Enable simplex solver plugin for VRM")
SET (USE_EXTERNAL_SUM "None" CACHE STRING "Enable external velocity solver, if any")
SET_PROPERTY(CACHE USE_EXTERNAL_SUM PROPERTY STRINGS "None" "onbody" "fmm2d")
SET (USE_HO_HYBRID "None" CACHE STRING "Select external hybrid solver, if any")
SET_PROPERTY(CACHE USE_HO_HYBRID PROPERTY STRINGS "None" "CXX" "Fortran")
SET (DO_PROFILE FALSE CACHE BOOL "Enable profiling via gprof")
SET (CMAKE_VERBOSE_MAKEFILE on)
SET (CMAKE_EXPORT_COMPILE_COMMANDS on)

# Output more info if debuging
IF (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
  ADD_DEFINITIONS(-DVERBOSE=true)
ELSE ()
  ADD_DEFINITIONS(-DVERBOSE=false)
ENDIF ()

IF (APPLE)
  SET (CMAKE_INSTALL_PREFIX /usr/local/share)
ENDIF ()

# Finds directory to pre-made sims
IF (${CMAKE_BUILD_TYPE} STREQUAL "Release")
  ADD_DEFINITIONS(-DEXAMPLES_DIR="${CMAKE_INSTALL_PREFIX}/examples")
ELSE ()
  ADD_DEFINITIONS(-DEXAMPLES_DIR="${CMAKE_CURRENT_SOURCE_DIR}/examples")
ENDIF ()


# Define compiler-specific options

IF (CMAKE_COMPILER_IS_GNUCXX)
  # support Core2 at a minimum, haswell ideally
  #SET (CMAKE_CXX_FLAGS "-Wall -Wformat -Wno-int-in-bool-context -std=c++17 -march=core2 -mtune=haswell")
  SET (CMAKE_CXX_FLAGS "-Wall -Wformat -std=c++17 -Wno-int-in-bool-context")
  SET (CMAKE_CXX_FLAGS_DEBUG "-O0 -g -ggdb3")
  SET (CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")
  SET (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -march=native -g -ggdb3")
  IF (DO_PROFILE)
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
  ENDIF ()
ELSEIF (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  SET (CMAKE_CXX_FLAGS "-Wall -Wformat -std=c++17")
  SET (CMAKE_CXX_FLAGS_DEBUG "-O0 -g -ggdb")
  SET (CMAKE_CXX_FLAGS_RELEASE "-O3")
  SET (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g -ggdb")
ELSEIF (MSVC)
  SET (CMAKE_CXX_FLAGS "/std:c++17 /EHsc /D_USE_MATH_DEFINES /DNOMINMAX")
  SET (CMAKE_CXX_FLAGS_DEBUG "/Zi")
  SET (CMAKE_CXX_FLAGS_RELEASE "/O2")
ENDIF ()


# Define platform-specific libraries and options

IF( WIN32 )
  # because clang is stupid
  IF (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
  ENDIF ()

  SET( GL_LIBS_HOME "C:/Lib" CACHE PATH "Location of support libraries" )
  SET( GLFW3_ROOT "${GL_LIBS_HOME}/glfw-3.3.bin.WIN32" CACHE STRING "GLFW3 root" )
  INCLUDE_DIRECTORIES( "${GLFW3_ROOT}/include" )
  LINK_DIRECTORIES   ( "${GLFW3_ROOT}/lib-vc2017" )
  SET( FRAMEWORK_LIBS glfw3.lib )

  SET( EIGEN_ROOT "${GL_LIBS_HOME}/eigen-3.3.7" CACHE STRING "Eigen root" )

  # OpenMP for multithreading
  IF( USE_OMP )
    FIND_PACKAGE( OpenMP REQUIRED )
    IF( OPENMP_FOUND )
      #message("OPENMP FOUND")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    ENDIF()
  ELSE()
    # suppress the warning
    #SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
    # need to find pthreads (for std::async?)
    #FIND_PACKAGE( Threads REQUIRED )
    #SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_THREAD_LIBS_INIT}")
    #IF (CMAKE_THREAD_LIBS_INIT)
    #ENDIF()
  ENDIF()
  
  # Vc for explicit vectorization
  IF( USE_VC )
    SET (CPREPROCDEFS ${CPREPROCDEFS} -DUSE_VC)
    SET (VC_ROOT "${GL_LIBS_HOME}/Vc" CACHE STRING "Vc root" )
    INCLUDE_DIRECTORIES( "${VC_ROOT}/include" )
    FIND_LIBRARY( VC_LIBS NAMES Vc PATHS "${VC_ROOT}/lib" )
  ELSE()
    SET (VC_LIBS "")
  ENDIF()
  
  SET( BASE_LIBS ${VC_LIBS} ${CMAKE_THREAD_LIBS_INIT} )
  SET( GUI_LIBS ${FRAMEWORK_LIBS} opengl32 )
  SET_PROPERTY(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT Omega2D)

ELSEIF( APPLE )
  # because clang is stupid
  IF (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
  ENDIF ()

  INCLUDE_DIRECTORIES ( /System/Library/Frameworks )
  FIND_LIBRARY(COCOA_LIBRARY Cocoa)
  FIND_LIBRARY(IOKIT_LIBRARY IOKit)
  FIND_LIBRARY(COREVIDEO_LIBRARY CoreVideo)
  FIND_LIBRARY(OpenGL_LIBRARY OpenGL )
  MARK_AS_ADVANCED (COCOA_LIBRARY
                    IOKIT_LIBRARY
                    OpenGL_LIBRARY)

  #SET( FRAMEWORK_LIBS glfw3 )
  FIND_LIBRARY( FRAMEWORK_LIBS NAMES glfw3 glfw libglfw PATHS /usr/local/lib)

  #INCLUDE_DIRECTORIES ( "/usr/local/include" )
  #INCLUDE_DIRECTORIES ( "/usr/local/include/eigen3" )
  #INCLUDE_DIRECTORIES ( "/usr/local/Cellar/eigen/3.3.4/include/eigen3" )
  SET( EIGEN_ROOT "/usr/local/include/eigen3" CACHE STRING "Eigen root" )

  IF( USE_OMP )
    FIND_PACKAGE( OpenMP )
    IF( OPENMP_FOUND )
      #message("OPENMP FOUND")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    ENDIF()
  ENDIF()

  # Vc for explicit vectorization
  IF( USE_VC )
    SET (CPREPROCDEFS ${CPREPROCDEFS} -DUSE_VC)
    SET (VC_ROOT "/opt/Vc" CACHE STRING "Vc root" )
    INCLUDE_DIRECTORIES( "${VC_ROOT}/include" )
    FIND_LIBRARY( VC_LIBS NAMES Vc PATHS "${VC_ROOT}/lib" )
  ELSE()
    SET (VC_LIBS "")
  ENDIF()

  SET( CMAKE_MACOSX_RPATH 1 )
  SET( BASE_LIBS ${VC_LIBS} ${CMAKE_THREAD_LIBS_INIT} )
  SET( GUI_LIBS ${FRAMEWORK_LIBS} ${COCOA_LIBRARY} ${IOKIT_LIBRARY} ${COREVIDEO_LIBRARY} ${OpenGL_LIBRARY} )

ELSEIF( UNIX )
  # because clang is stupid
  IF (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++")
  ENDIF ()

  # glfw finds opengl
  # this works for CentOS 7
  #FIND_PACKAGE( glfw3 3.0 REQUIRED )
  # this works for Fedora 21
  #FIND_PACKAGE( glfw 3.0 REQUIRED )
  # this seems to be more general
  FIND_LIBRARY( FRAMEWORK_LIBS NAMES glfw3 glfw libglfw PATHS /usr/local/lib )

  # tell imgui to use glad
  SET (CPREPROCDEFS ${CPREPROCDEFS} -DIMGUI_IMPL_OPENGL_LOADER_GLAD)

  # but glad needs dl, which cmake doesn't think glfw needs
  SET( FRAMEWORK_LIBS glfw dl )

  SET( EIGEN_ROOT "/usr/include/eigen3" CACHE STRING "Eigen root" )

  # OpenMP for multithreading
  IF( USE_OMP )
    FIND_PACKAGE( OpenMP REQUIRED )
    IF( OPENMP_FOUND )
      #message("OPENMP FOUND")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    ENDIF()
  ELSE()
    # suppress the warning
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
    # need to find pthreads (for std::async?)
    FIND_PACKAGE( Threads REQUIRED )
    #SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_THREAD_LIBS_INIT}")
    #IF (CMAKE_THREAD_LIBS_INIT)
    #ENDIF()
  ENDIF()

  # Vc for explicit vectorization
  IF( USE_VC )
    SET (CPREPROCDEFS ${CPREPROCDEFS} -DUSE_VC)
    SET (VC_ROOT "/opt/Vc" CACHE STRING "Vc root" )
    INCLUDE_DIRECTORIES( "${VC_ROOT}/include" )
    FIND_LIBRARY( VC_LIBS NAMES Vc PATHS "${VC_ROOT}/lib" )
  ELSE()
    SET (VC_LIBS "")
  ENDIF()

  SET( BASE_LIBS ${VC_LIBS} ${CMAKE_THREAD_LIBS_INIT} )
  SET( GUI_LIBS ${FRAMEWORK_LIBS} ${GLFW_LIBRARIES} )
ENDIF()


# Enable plugins

# experimental adaptive VRM
IF( USE_PLUGIN_AVRM )
  SET (CPREPROCDEFS ${CPREPROCDEFS} -DPLUGIN_AVRM)
  INCLUDE_DIRECTORIES( "plugin-o2d-avrm" )
  # add VRMadaptive.h to HEADERS?
ENDIF()

# simplex solver
IF( USE_PLUGIN_SIMPLEX )
  SET (CPREPROCDEFS ${CPREPROCDEFS} -DPLUGIN_SIMPLEX)
  INCLUDE_DIRECTORIES( "plugin-o2d-simplex" )
  # add simplex.h to HEADERS?
ENDIF()

# external fast summation routine
IF (${USE_EXTERNAL_SUM} STREQUAL "onbody")
  ADD_SUBDIRECTORY( extern/onbody )
  INCLUDE_DIRECTORIES( "extern/onbody/src" )
  SET (EXTERNAL_LIBS ${EXTERNAL_LIBS} bh2dvort)
  SET (CPREPROCDEFS ${CPREPROCDEFS} -DEXTERNAL_VEL_SOLVE)
ELSEIF( ${USE_EXTERNAL_SUM} STREQUAL "fmm2d" )
  SET (CPREPROCDEFS ${CPREPROCDEFS} -DEXTERNAL_VEL_SOLVE)
  SET (FASTSUM_LIBS "libfastsumm.a" CACHE STRING "Static library exposing velocity summation function call" )
  SET (EXTERNAL_LIBS ${FASTSUM_LIBS} gfortran)
ELSE()
  SET (FASTSUM_LIBS "")
ENDIF()

# Enable submodules
ADD_SUBDIRECTORY( extern/gmsh-reader )
INCLUDE_DIRECTORIES( "extern/gmsh-reader/src" )
SET( EXTERNAL_LIBS ${EXTERNAL_LIBS} gmsh-reader-static )

IF (${USE_HO_HYBRID} STREQUAL "Fortran")
  ADD_SUBDIRECTORY( extern/HO-Fortran )
  INCLUDE_DIRECTORIES( "extern/HO-Fortran/src" )
  SET (EXTERNAL_LIBS ${EXTERNAL_LIBS} 2DVortTrans ${APLLES_LIBS} ${XML_LIBS} ${LAPACK_LIBS} gfortran )
  SET (CPREPROCDEFS ${CPREPROCDEFS} -DHOFORTRAN)
ELSEIF( ${USE_HO_HYBRID} STREQUAL "CXX" )
  ADD_SUBDIRECTORY( extern/HO-CXX )
  INCLUDE_DIRECTORIES ("extern/HO-CXX/src" )
  #INCLUDE_DIRECTORIES ( "${EIGEN_ROOT}" )
  INCLUDE_DIRECTORIES ( "${AMGCL_ROOT}" )
  SET (EXTERNAL_LIBS ${EXTERNAL_LIBS} HO-CXX)
  SET (CPREPROCDEFS ${CPREPROCDEFS} -DHOCXX)
ENDIF()

# Define files and targets

ADD_DEFINITIONS (${CPREPROCDEFS})

INCLUDE_DIRECTORIES( "${EIGEN_ROOT}" )
INCLUDE_DIRECTORIES( "lib" )
# For some of the libraries who need glad
INCLUDE_DIRECTORIES( "lib/glad" )
# For some of the internal libraries who need nanoflann
INCLUDE_DIRECTORIES( "lib/nanoflann" )
# For some of the internal libraries who need nnls 
INCLUDE_DIRECTORIES( "lib/eigen-nnls" )
INCLUDE_DIRECTORIES( "src" )

IF (CMAKE_COMPILER_IS_GNUCXX)
  SET_SOURCE_FILES_PROPERTIES("lib/imgui/imgui_draw.cpp" PROPERTIES COMPILE_FLAGS -Wno-maybe-uninitialized)
ENDIF()

FILE(GLOB IMGUISOURCES "lib/imgui/*.cpp" )
SET(SOURCES "src/Simulation.cpp"
            "src/FlowFeature.cpp"
            "src/CollectionHelper.cpp"
            "src/BoundaryFeature.cpp"
            "src/MeasureFeature.cpp"
            "src/Body.cpp"
            "src/SegmentHelper.cpp"
            "src/RenderParams.cpp"
            "src/JsonHelper.cpp"
            "src/StatusFile.cpp"
            "lib/tinyxml2/tinyxml2.cpp"
            "lib/tinyexpr/tinyexpr.c" )
SET(GUI_SOURCES "lib/glad/glad.c"
                "src/ShaderHelper.cpp"
                "src/FeatureDraw.cpp"
                "lib/miniz/miniz.c"
                "lib/miniz/FrameBufferToImage.cpp" )

# create a binary for the GUI version
IF( BUILD_GUI )
  ADD_EXECUTABLE( ${PROJECT_NAME} ${IMGUISOURCES} ${SOURCES} ${GUI_SOURCES} "src/main_gui.cpp" )
  SET_TARGET_PROPERTIES( ${PROJECT_NAME} PROPERTIES OUTPUT_NAME "${PROJECT_NAME}.bin" )
  TARGET_LINK_LIBRARIES( ${PROJECT_NAME} ${BASE_LIBS} ${GUI_LIBS} ${EXTERNAL_LIBS} )
  TARGET_COMPILE_DEFINITIONS( ${PROJECT_NAME} PRIVATE "-DUSE_GL" )
  TARGET_COMPILE_DEFINITIONS( ${PROJECT_NAME} PRIVATE "-DUSE_IMGUI" )
  INSTALL( TARGETS ${PROJECT_NAME} DESTINATION bin )
ENDIF()

# create a binary for the command-line version
IF( BUILD_BATCH )
  ADD_EXECUTABLE( "${PROJECT_NAME}batch" ${SOURCES} "src/main_batch.cpp" )
  SET_TARGET_PROPERTIES( "${PROJECT_NAME}batch" PROPERTIES OUTPUT_NAME "${PROJECT_NAME}batch.bin" )
  TARGET_LINK_LIBRARIES( "${PROJECT_NAME}batch" ${BASE_LIBS} ${EXTERNAL_LIBS} )
  INSTALL( TARGETS "${PROJECT_NAME}batch" DESTINATION bin )
ENDIF()

INSTALL( DIRECTORY examples/ DESTINATION examples )
INSTALL( FILES LICENSE DESTINATION LICENSE )

