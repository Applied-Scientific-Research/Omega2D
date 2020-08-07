//
// FrameBufferToImage.cpp - header for calling glReadPixels and stb_image_write
//                          parts from https://github.com/gileoo/Imgui-IGS-Snippets
//

#include "FrameBufferToImage.h"

// header-only png writing
/*#define STB_IMAGE_WRITE_IMPLEMENTATION
#define STBI_MSC_SECURE_CRT
#include "stb_image_write.h"*/

#define MINIZ_NO_STDIO
#define MINIZ_NO_TIME
#define MINIZ_NO_ZLIB_APIS
#include "miniz.h"
#include <stdio.h>
#include <vector>

void saveFramePNG(std::string file_name) {

  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT, viewport);

  // assure alignment of 4, png size needs to by dividable by 4
  int width  = (viewport[2]/4) * 4;
  int height = (viewport[3]/4) * 4;
  int nPixels = width * height;

  // GL_RGB implies "3" in two other places

  GLubyte* buffer = new GLubyte[nPixels*3];
	
  glReadBuffer(GL_BACK);
  glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, buffer);

  //stbi_flip_vertically_on_write(1);
  //int stbi_write_png(char const *filename, int w, int h, int comp, const void *data, int stride_in_bytes);
  //stbi_write_png(file_name.c_str(), width, height, 3, buffer, 3*width);

  size_t png_data_size = 0;
  void *pPNG_data = tdefl_write_image_to_png_file_in_memory_ex(buffer, width, height, 3,
                                                               &png_data_size, MZ_BEST_COMPRESSION, MZ_TRUE);
  if (pPNG_data) {
    FILE *pFile = fopen(file_name.c_str(), "wb");
    fwrite(pPNG_data, 1, png_data_size, pFile);
    fclose(pFile);
    printf("Wrote %s\n", file_name.c_str());
  } else {
    fprintf(stderr, "tdefl_write_image_to_png_file_in_memory_ex() failed!\n");
  }

  mz_free(pPNG_data);
  std::free(buffer);
}
