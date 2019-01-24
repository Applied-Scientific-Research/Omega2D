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
#include <cstdio>	// for FILE
#include <sstream>	// for stringstream
#include <iomanip>	// for setfill, setw

//
// write point data to a .vtu file
//
// should we be using PolyData or UnstructuredGrid for particles?
//
template <class S>
void write_vtu_points(Points<S> const& pts) {

  // generate file name
  static int frameno = 0;
  std::stringstream vtkfn;
  vtkfn << "part_" << std::setfill('0') << std::setw(5) << frameno << ".vtu";

  // prepare file pointer and printer
  std::FILE* fp = std::fopen(vtkfn.str().c_str(), "w");
  tinyxml2::XMLPrinter printer( fp );

  // write <?xml version="1.0"?>
  printer.PushHeader(false, true);

  printer.OpenElement( "VTKFile" );
  printer.PushAttribute( "type", "UnstructuredGrid" );
  //printer.PushAttribute( "type", "PolyData" );
  printer.PushAttribute( "version", "0.1" );
  printer.PushAttribute( "byte_order", "LittleEndian" );

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
  printer.PushAttribute( "format", "ascii" );
  // end the element header and write the data
  printer.PushText( " " );
  //printer.PushText( "" );
  const std::array<Vector<S>,Dimensions>& x = pts.get_pos();
  for (size_t i=0; i<pts.get_n(); ++i) {
    std::fprintf(fp, "%g %g 0.0 ", x[0][i], x[1][i]);
  }
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
  printer.PushAttribute( "Scalars", "circulation" );

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
  printer.CloseElement();	// DataArray

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

  printer.CloseElement();	// Piece
  printer.CloseElement();	// PolyData or UnstructuredGrid
  printer.CloseElement();	// VTKFile

  std::fclose(fp);

  // increment and return
  frameno++;
}
