//
// FrameBufferToImage.h - header for calling glReadPixels and stb_image_write
//                        from https://github.com/gileoo/Imgui-IGS-Snippets
//

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


void saveFramePNG(std::string file_name);

