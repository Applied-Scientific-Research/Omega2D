#
# dummy_solver
# (c)2020 Applied Scientific Research, Inc.
#
CMAKE_MINIMUM_REQUIRED( VERSION 2.8 )
PROJECT( dummy_solver )
ENABLE_LANGUAGE (CXX)

# User-serviceable variables here
SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo")
IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE "Release")
ENDIF ()
SET (CMAKE_VERBOSE_MAKEFILE on)
SET (CMAKE_EXPORT_COMPILE_COMMANDS on)

# Define compiler-specific options
IF (CMAKE_COMPILER_IS_GNUCXX)
  SET (CMAKE_CXX_FLAGS "-Wall -Wformat -std=c++17 -Wno-int-in-bool-context")
  SET (CMAKE_CXX_FLAGS_DEBUG "-O0 -g -ggdb3")
  SET (CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")
  SET (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -march=native -g -ggdb3")
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

# Define targets
ADD_LIBRARY ( "${PROJECT_NAME}-static" STATIC "dummysolver.cpp" )
ADD_LIBRARY ( ${PROJECT_NAME} SHARED "dummysolver.cpp" )

#ADD_EXECUTABLE( "${PROJECT_NAME}.bin" "main.cpp" )
#TARGET_LINK_LIBRARIES( "${PROJECT_NAME}.bin" LINK_PUBLIC "${PROJECT_NAME}-static" )
