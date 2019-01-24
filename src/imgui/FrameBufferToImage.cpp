#include <vector>
#include "FrameBufferToImage.hpp"

using namespace std;

void writePngFile(const char *filename, int width, int height, GLubyte* pixels )
{
	FILE *fp = fopen(filename, "wb");
	if (!fp) 
		abort();

	png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	if (!png) 
		abort();

	png_infop info = png_create_info_struct(png);
	if (!info) abort();

	if (setjmp(png_jmpbuf(png))) abort();

	png_init_io(png, fp);

	// Output is 8bit depth, RGB format.
	png_set_IHDR(
		png,
		info,
		width, height,
		8,
		PNG_COLOR_TYPE_RGB,
		PNG_INTERLACE_NONE,
		PNG_COMPRESSION_TYPE_DEFAULT,
		PNG_FILTER_TYPE_BASE
	);

	png_colorp palette = (png_colorp)png_malloc(png, PNG_MAX_PALETTE_LENGTH * sizeof(png_color));
	if (!palette) 
	{
		fclose(fp);
		png_destroy_write_struct(&png, &info);
		return;
	}
	png_set_PLTE(png, info, palette, PNG_MAX_PALETTE_LENGTH);
	png_write_info(png, info);
	png_set_packing(png);

	png_bytepp rows = (png_bytepp)png_malloc(png, height * sizeof(png_bytep));
	for (int i = 0; i < height; ++i)
		rows[i] = (png_bytep)(pixels + (height - i - 1) * width * 3);

	png_write_image(png, rows);
	
	png_write_end(png, info);
	png_free(png, palette);
	png_destroy_write_struct(&png, &info);
	   
	free(rows);
	fclose(fp);
}

void saveFramePNG( string file_name )
{
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);

	// assure alignment of 4, png size needs to by dividable by 4
	int width  = (viewport[2]/4) * 4;
	int height = (viewport[3]/4) * 4;
	int nPixels = width * height;

	GLubyte* buffer = new GLubyte[nPixels*3];
	
	glReadBuffer(GL_BACK);
	glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, buffer);

	writePngFile(file_name.c_str(), width, height, buffer);

	free(buffer);
}

