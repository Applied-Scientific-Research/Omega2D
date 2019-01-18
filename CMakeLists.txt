#
# Omega2D
# (c)2017-2018 Applied Scientific Research, Inc.
#
CMAKE_MINIMUM_REQUIRED( VERSION 2.8 )
PROJECT( Omega2D )

SET (CMAKE_CONFIGURATION_TYPES "Debug" "Release")
#if(NOT CMAKE_BUILD_TYPE)
#    set(CMAKE_BUILD_TYPE Release)
#endif(NOT CMAKE_BUILD_TYPE)
SET (CMAKE_BUILD_TYPE "Release" CACHE STRING "Select which configuration to build" )
SET (USE_OMP FALSE CACHE BOOL "Use OpenMP multithreading")
SET (USE_VC FALSE CACHE BOOL "Use Vc for vector arithmetic")
SET (CMAKE_VERBOSE_MAKEFILE on)
SET (CMAKE_EXPORT_COMPILE_COMMANDS on)

IF (CMAKE_COMPILER_IS_GNUCXX)
  SET (CMAKE_CXX_FLAGS "-Wall -Wformat -std=c++17 -Wno-int-in-bool-context")
  SET (CMAKE_CXX_FLAGS_DEBUG "-O0 -g -ggdb")
  SET (CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")
ELSEIF (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  SET (CMAKE_CXX_FLAGS "-Wall -Wformat -std=c++17 -stdlib=libc++")
  SET (CMAKE_CXX_FLAGS_DEBUG "-O0 -g -ggdb")
  SET (CMAKE_CXX_FLAGS_RELEASE "-O3")
ELSEIF (MSVC)
  SET (CMAKE_CXX_FLAGS "/Ox")
ENDIF ()

IF( WIN32 )
  SET( GL_LIBS_HOME "C:/lib" CACHE PATH "Location of support libraries" )
  SET( GLFW3_ROOT "${GL_LIBS_HOME}/glfw-3.1.2.bin.WIN32" CACHE STRING "GLFW3 root" )
  INCLUDE_DIRECTORIES( "${GLFW3_ROOT}/include" )
  LINK_DIRECTORIES   ( "${GLFW3_ROOT}/lib-vc2015" )
  SET( FRAMEWORK_LIBS glfw3.lib )
  SET( PLATFORM_LIBS ${FRAMEWORK_LIBS} opengl32 )

ELSEIF( APPLE )
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

  INCLUDE_DIRECTORIES ( "/usr/local/include" )
  INCLUDE_DIRECTORIES ( "/usr/local/include/eigen3" )
  #INCLUDE_DIRECTORIES ( "/usr/local/Cellar/eigen/3.3.4/include/eigen3" )

  SET(PLATFORM_LIBS ${FRAMEWORK_LIBS} ${COCOA_LIBRARY} ${IOKIT_LIBRARY} ${COREVIDEO_LIBRARY} ${OpenGL_LIBRARY})

  IF( USE_OMP )
    FIND_PACKAGE( OpenMP )
    IF( OPENMP_FOUND )
      #message("OPENMP FOUND")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    ENDIF()
  ENDIF()

ELSEIF( UNIX )
  # glfw finds opengl
  # this works for CentOS 7
  #FIND_PACKAGE( glfw3 3.0 REQUIRED )
  # this works for Fedora 21
  #FIND_PACKAGE( glfw 3.0 REQUIRED )
  # this seems to be more general
  FIND_LIBRARY( FRAMEWORK_LIBS NAMES glfw3 glfw libglfw PATHS /usr/local/lib )

  # but glad needs dl, which cmake doesn't think glfw needs
  SET( FRAMEWORK_LIBS glfw dl )

  INCLUDE_DIRECTORIES ( "/usr/include/eigen3" )

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
    INCLUDE_DIRECTORIES( "/opt/Vc/include" )
    FIND_LIBRARY( VC_LIBS NAMES Vc PATHS /opt/Vc/lib )
  ELSE()
    SET (VC_LIBS "")
  ENDIF()

  SET( PLATFORM_LIBS ${FRAMEWORK_LIBS} ${GLFW_LIBRARIES} ${VC_LIBS} ${CMAKE_THREAD_LIBS_INIT})
ENDIF()

ADD_DEFINITIONS (${CPREPROCDEFS})
set(GLADSRC "src/glad/glad.c" )
INCLUDE_DIRECTORIES( "src/glad" )
SET(SOURCES "src/Simulation.cpp"
            "src/FlowFeature.cpp"
            "src/CollectionHelper.cpp"
            "src/BoundaryFeature.cpp"
            "src/MeasureFeature.cpp"
            "src/main.cpp")
file(GLOB IMGUISOURCES "src/imgui/*.cpp" )

# create a binary for the GUI version
ADD_EXECUTABLE( ${PROJECT_NAME} ${GLADSRC} ${IMGUISOURCES} ${SOURCES} "src/ShaderHelper.cpp" )
SET_TARGET_PROPERTIES( ${PROJECT_NAME} PROPERTIES OUTPUT_NAME "${PROJECT_NAME}.bin" )
TARGET_LINK_LIBRARIES( ${PROJECT_NAME} ${PLATFORM_LIBS} )
TARGET_COMPILE_DEFINITIONS( ${PROJECT_NAME} PRIVATE "-DUSE_GL" )

