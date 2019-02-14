/*
 * VtkXmlHelper.h - Write an XML-format VTK data file using TinyXML2
 *
 * (c)2019 Applied Scientific Research, Inc.
 *         Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Omega2D.h"
#include "VectorHelper.h"
#include "Points.h"

#include "tinyxml2.h"

#include <vector>
#include <cstdint>	// for uint32_t
#include <cstdio>	// for FILE
#include <sstream>	// for stringstream
#include <iomanip>	// for setfill, setw

//
// compress a byte stream
//


//
// convert an array or interleave arrays to generate a base64 bytestream
//
// why would you ever want to use base64 for floats and such? so wasteful.
//
template <class S>
const char* convert_to_base64 (Vector<S> const & v, bool compress = false) {
  const char* retval = "hi";

  if (compress) {
  }

  return retval;
}

//const char* interleave_to_base64 (Vector<S>


//
// write point data to a .vtu file
//
// should we be using PolyData or UnstructuredGrid for particles?
//
// note bad documetation, somewhat corrected here:
// https://mathema.tician.de/what-they-dont-tell-you-about-vtk-xml-binary-formats/
//
template <class S>
void write_vtu_points(Points<S> const& pts, const size_t file_idx, const size_t frameno) {

  bool has_radii = true;
  bool has_strengths = true;
  std::string prefix = "part_";
  if (pts.is_inert()) {
    has_strengths = false;
    has_radii = false;
    prefix = "fldpt_";
  }

  // generate file name
  std::stringstream vtkfn;
  vtkfn << prefix << std::setfill('0') << std::setw(2) << file_idx << "_" << std::setw(5) << frameno << ".vtu";

  // prepare file pointer and printer
  std::FILE* fp = std::fopen(vtkfn.str().c_str(), "wb");
  tinyxml2::XMLPrinter printer( fp );

  // write <?xml version="1.0"?>
  printer.PushHeader(false, true);

  printer.OpenElement( "VTKFile" );
  printer.PushAttribute( "type", "UnstructuredGrid" );
  //printer.PushAttribute( "type", "PolyData" );
  printer.PushAttribute( "version", "0.1" );
  printer.PushAttribute( "byte_order", "LittleEndian" );
  printer.PushAttribute( "header_type", "UInt32" );

  // push comment with sim time?

  printer.OpenElement( "UnstructuredGrid" );
  //printer.OpenElement( "PolyData" );
  printer.OpenElement( "Piece" );
  printer.PushAttribute( "NumberOfPoints", std::to_string(pts.get_n()).c_str() );
  printer.PushAttribute( "NumberOfCells", std::to_string(pts.get_n()).c_str() );

  printer.OpenElement( "Points" );
  printer.OpenElement( "DataArray" );
  printer.PushAttribute( "NumberOfComponents", "3" );
  printer.PushAttribute( "Name", "position" );
  printer.PushAttribute( "type", "Float32" );
  const std::array<Vector<S>,Dimensions>& x = pts.get_pos();
  // write ascii version
  printer.PushAttribute( "format", "ascii" );
  printer.PushText( " " );
  for (size_t i=0; i<pts.get_n(); ++i) { std::fprintf(fp, "%g %g 0.0 ", x[0][i], x[1][i]); }
  // write binary version
  //printer.PushAttribute( "format", "binary" );
  //printer.PushText( "\n" );
  //printer.PrintSpace(5);
  //std::fwrite(v.data(), sizeof v[0], v.size(), f1);
  //const float zero = 0.0;
  //for (size_t i=0; i<pts.get_n(); ++i) {
  //  std::fwrite(&x[0][i], sizeof(S), 1, fp);
  //  std::fwrite(&x[1][i], sizeof(S), 1, fp);
  //  std::fwrite(&zero, sizeof(S), 1, fp);
  //}
  //printer.PushText( "\n" );
  //printer.PrintSpace(5);
  // close it out
  printer.CloseElement();	// DataArray
  printer.CloseElement();	// Points

  printer.OpenElement( "Cells" );

  printer.OpenElement( "DataArray" );
  printer.PushAttribute( "Name", "connectivity" );
  printer.PushAttribute( "type", "UInt32" );
  printer.PushAttribute( "format", "ascii" );
  printer.PushText( " " );
  for (size_t i=0; i<pts.get_n(); ++i) { std::fprintf(fp, "%d ", (int)i); }
  printer.CloseElement();	// DataArray

  printer.OpenElement( "DataArray" );
  printer.PushAttribute( "Name", "offsets" );
  printer.PushAttribute( "type", "UInt32" );
  printer.PushAttribute( "format", "ascii" );
  printer.PushText( " " );
  for (size_t i=0; i<pts.get_n(); ++i) { std::fprintf(fp, "%d ", (int)i); }
  printer.CloseElement();	// DataArray

  printer.OpenElement( "DataArray" );
  printer.PushAttribute( "Name", "types" );
  printer.PushAttribute( "type", "UInt8" );
  printer.PushAttribute( "format", "ascii" );
  printer.PushText( " " );
  for (size_t i=0; i<pts.get_n(); ++i) { std::fprintf(fp, "%d ", (int)1); }
  printer.CloseElement();	// DataArray

  printer.CloseElement();	// Cells

  printer.OpenElement( "PointData" );
  printer.PushAttribute( "Vectors", "velocity" );
  std::string scalar_list;
  if (has_strengths) scalar_list.append("circulation,");
  if (has_radii) scalar_list.append("radius,");
  if (scalar_list.size()>1) {
    scalar_list.pop_back();
    printer.PushAttribute( "Scalars", scalar_list.c_str() );
  }

  if (has_strengths) {
    printer.OpenElement( "DataArray" );
    printer.PushAttribute( "Name", "circulation" );
    printer.PushAttribute( "type", "Float32" );
    printer.PushAttribute( "format", "ascii" );
    // end the element header and write the data
    printer.PushText( " " );
    const Vector<S>& s = pts.get_str();
    for (size_t i=0; i<pts.get_n(); ++i) {
      std::fprintf(fp, "%g ", s[i]);
    }
    // write binary version
    //printer.PushAttribute( "format", "appended" );
    //printer.PushAttribute( "offset", "0" );
    printer.CloseElement();	// DataArray
  }

  if (has_radii) {
    printer.OpenElement( "DataArray" );
    printer.PushAttribute( "Name", "radius" );
    printer.PushAttribute( "type", "Float32" );
    printer.PushAttribute( "format", "ascii" );
    // end the element header and write the data
    printer.PushText( " " );
    const Vector<S>& r = pts.get_rad();
    for (size_t i=0; i<pts.get_n(); ++i) {
      std::fprintf(fp, "%g ", r[i]);
    }
    printer.CloseElement();	// DataArray
  }

  printer.OpenElement( "DataArray" );
  printer.PushAttribute( "NumberOfComponents", "3" );
  printer.PushAttribute( "Name", "velocity" );
  printer.PushAttribute( "type", "Float32" );
  printer.PushAttribute( "format", "ascii" );
  // end the element header and write the data
  printer.PushText( " " );
  const std::array<Vector<S>,Dimensions>& u = pts.get_vel();
  for (size_t i=0; i<pts.get_n(); ++i) {
    std::fprintf(fp, "%g %g 0.0 ", u[0][i], u[1][i]);
  }
  printer.CloseElement();	// DataArray

  printer.CloseElement();	// PointData

  //printer.OpenElement( "CellData" );
  //printer.CloseElement();	// CellData

  // here's the problem: ParaView's VTK/XML reader is not able to read "raw" bytestreams
  // it parses each character and inevitably sees a > or <, so complains about mismatched
  //   tags, and doesn't seem to point to the right place
  // everyone who complains about raw data in the AppendedData section makes the mistake
  //   of forgetting the 4-byte length - I've got that
  // it's VTK's fault here
  // Jesus, for how inefficient saving 2D particle data is, it's compounded by having to 
  //   do it all in ascii!!! VTK just sucks. That's really what my 6 hours of wasted work
  //   has taught me. It just sucks. Don't use it. Not that anything else is better.
/*
  printer.OpenElement( "AppendedData" );
  printer.PushAttribute( "encoding", "raw" );
  printer.PushText( " " );
  printer.PushText( "_" );
  uint32_t arry_len = s.size() * sizeof(s[0]);
  char* ptr = (char*)(&arry_len);
  std::cout << "  writing " << arry_len << " bytes to appended data" << std::endl;
  //std::fwrite(&arry_len, sizeof(uint32_t), 1, fp);
  std::fwrite(ptr, sizeof(uint32_t), 1, fp);
  std::fwrite(s.data(), sizeof(s[0]), s.size(), fp);
  printer.PushText( " " );
  printer.CloseElement();	// AppendedData
*/

  printer.CloseElement();	// Piece
  printer.CloseElement();	// PolyData or UnstructuredGrid
  printer.CloseElement();	// VTKFile

  std::fclose(fp);

  std::cout << "Wrote particle data to " << vtkfn.str() << std::endl;
}
