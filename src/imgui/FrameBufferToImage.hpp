#pragma once

#include <string>

//#include "GL/glew.h"
#ifdef _WIN32
  // for glad
  #define APIENTRY __stdcall
  // for C++11 stuff
  #include <ciso646>
#endif
#include "glad.h"

//http://www.libpng.org/pub/png/libpng.html
#include "png.h"

void writePngFile(const char *filename, int width, int height, GLubyte* pixels );
void saveFramePNG( std::string file_name );

