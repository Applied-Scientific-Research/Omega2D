/*
 * VtkXmlHelper.h - Write an XML-format VTK data file using TinyXML2
 *
 * (c)2019 Applied Scientific Research, Inc.
 *         Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Omega2D.h"
#include "VectorHelper.h"
#include "Collection.h"
#include "Points.h"
#include "Surfaces.h"

#include "tinyxml2.h"
#include "cppcodec/base64_rfc4648.hpp"

#include <vector>
#include <cstdint>	// for uint32_t
#include <cstdio>	// for FILE
#include <string>	// for string
#include <sstream>	// for stringstream
#include <iomanip>	// for setfill, setw


//
// compress a byte stream - TO DO
//


//
// write vector to the vtk file
//
// why would you ever want to use base64 for floats and such? so wasteful.
//
template <class S>
void write_DataArray (tinyxml2::XMLPrinter& _p,
                      Vector<S> const & _data,
                      bool _compress = false,
                      bool _asbase64 = false) {

  using base64 = cppcodec::base64_rfc4648;

  if (_compress) {
    // not yet implemented
  }

  if (_asbase64) {
    _p.PushAttribute( "format", "binary" );

    std::string encoded = base64::encode((const char*)_data.data(), (size_t)(sizeof(_data[0])*_data.size()));

    // encode and write the UInt32 length indicator - the length of the base64-encoded blob
    uint32_t encoded_length = encoded.size();
    std::string header = base64::encode((const char*)(&encoded_length), (size_t)(sizeof(uint32_t)));
    //std::cout << "base64 length of header: " << header.size() << std::endl;
    //std::cout << "base64 length of data: " << encoded.size() << std::endl;

    _p.PushText( " " );
    _p.PushText( header.c_str() );
    _p.PushText( encoded.c_str() );
    _p.PushText( " " );

    if (false) {
      std::vector<uint8_t> decoded = base64::decode(encoded);
      S* decoded_floats = reinterpret_cast<S*>(decoded.data());
  
      std::cout << "Encoding an array of " << _data.size() << " floats to: " << encoded << std::endl;
      std::cout << "Decoding back into array of floats, starting with: " << decoded_floats[0] << " " << decoded_floats[1] << std::endl;
      //std::cout << "Decoding back into array of " << decoded.size() << " floats, starting with: " << decoded[0] << std::endl;
    }

  } else {
    _p.PushAttribute( "format", "ascii" );

    // write the data
    _p.PushText( " " );
    for (size_t i=0; i<_data.size(); ++i) {
      _p.PushText( _data[i] );
      _p.PushText( " " );
    }
  }

  return;
}


//
// interleave arrays and write to the vtk file
//
template <class S>
void write_DataArray (tinyxml2::XMLPrinter& _p,
                      std::array<Vector<S>,2> const & _data,
                      bool _compress = false,
                      bool _asbase64 = false) {

  // interleave the two vectors into a new one
  Vector<S> newvec;
  newvec.resize(3 * _data[0].size());
  for (size_t i=0; i<_data[0].size(); ++i) {
    newvec[3*i+0] = _data[0][i];
    newvec[3*i+1] = _data[1][i];
    newvec[3*i+2] = 0.0;
  }

  // pass it on to write
  write_DataArray<S>(_p, newvec, _compress, _asbase64);
}


//
// interleave arrays and write to the vtk file
//
template <class S>
void write_DataArray (tinyxml2::XMLPrinter& _p,
                      std::array<Vector<S>,3> const & _data,
                      bool _compress = false,
                      bool _asbase64 = false) {

  // interleave the two vectors into a new one
  Vector<S> newvec;
  newvec.resize(3 * _data[0].size());
  for (size_t i=0; i<_data[0].size(); ++i) {
    newvec[3*i+0] = _data[0][i];
    newvec[3*i+1] = _data[1][i];
    newvec[3*i+2] = _data[2][i];
  }

  // pass it on to write
  write_DataArray<S>(_p, newvec, _compress, _asbase64);
}



//
// write point data to a .vtu file
//
// should we be using PolyData or UnstructuredGrid for particles?
//
// note bad documentation, somewhat corrected here:
// https://mathema.tician.de/what-they-dont-tell-you-about-vtk-xml-binary-formats/
//
// see the full vtk spec here:
// https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
//
// time can be written to a vtk file:
// https://gitlab.kitware.com/vtk/vtk/commit/6e0acf3b773f120c3b8319c4078a4eac9ed31ce1
//
// <PolyData>
//      <FieldData>
//        <DataArray type="Float64" Name="TimeValue" NumberOfTuples="1">1.24
//        </DataArray>
//      </FieldData>

//
template <class S>
std::string write_vtu_points(Points<S> const& pts, const size_t file_idx,
                             const size_t frameno, const double time) {

  assert(pts.get_n() > 0 && "Inside write_vtu_points with no points");

  const bool compress = false;
  const bool asbase64 = true;

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
  // note this is still unsigned even though all indices later are signed! thanks, Obama.
  printer.PushAttribute( "header_type", "UInt32" );

  // push comment with sim time?

  // must choose one of these two formats
  printer.OpenElement( "UnstructuredGrid" );
  //printer.OpenElement( "PolyData" );

  // include simulation time here
  printer.OpenElement( "FieldData" );
  printer.OpenElement( "DataArray" );
  printer.PushAttribute( "type", "Float64" );
  printer.PushAttribute( "Name", "TimeValue" );
  printer.PushAttribute( "NumberOfTuples", "1" );
  {
    Vector<double> time_vec = {time};
    write_DataArray (printer, time_vec, false, false);
  }
  printer.CloseElement();	// DataArray
  printer.CloseElement();	// FieldData

  printer.OpenElement( "Piece" );
  printer.PushAttribute( "NumberOfPoints", std::to_string(pts.get_n()).c_str() );
  printer.PushAttribute( "NumberOfCells", std::to_string(pts.get_n()).c_str() );

  printer.OpenElement( "Points" );
  printer.OpenElement( "DataArray" );
  printer.PushAttribute( "NumberOfComponents", "3" );
  printer.PushAttribute( "Name", "position" );
  printer.PushAttribute( "type", "Float32" );
  write_DataArray (printer, pts.get_pos(), compress, asbase64);
  printer.CloseElement();	// DataArray
  printer.CloseElement();	// Points

  printer.OpenElement( "Cells" );

  // https://discourse.paraview.org/t/cannot-open-vtu-files-with-paraview-5-8/3759
  // apparently the Vtk format documents indicate that connectivities and offsets
  //   must be in Int32, not UIntAnything. Okay...
  printer.OpenElement( "DataArray" );
  printer.PushAttribute( "Name", "connectivity" );
  printer.PushAttribute( "type", "Int32" );
  {
    Vector<int32_t> v(pts.get_n());
    std::iota(v.begin(), v.end(), 0);
    write_DataArray (printer, v, compress, asbase64);
  }
  printer.CloseElement();	// DataArray

  printer.OpenElement( "DataArray" );
  printer.PushAttribute( "Name", "offsets" );
  printer.PushAttribute( "type", "Int32" );
  {
    Vector<int32_t> v(pts.get_n());
    std::iota(v.begin(), v.end(), 1);
    write_DataArray (printer, v, compress, asbase64);
  }
  printer.CloseElement();	// DataArray

  // what about these? Int32 as well? hello? is there anybody there?
  printer.OpenElement( "DataArray" );
  printer.PushAttribute( "Name", "types" );
  printer.PushAttribute( "type", "UInt8" );
  {
    Vector<uint8_t> v(pts.get_n());
    std::fill(v.begin(), v.end(), 1);
    write_DataArray (printer, v, compress, asbase64);
  }
  printer.CloseElement();	// DataArray

  printer.CloseElement();	// Cells


  printer.OpenElement( "PointData" );
  std::string vector_list;
  std::string scalar_list;

  vector_list.append("velocity,");
  if (has_strengths) scalar_list.append("circulation,");
  if (has_radii) scalar_list.append("radius,");

  if (vector_list.size()>1) {
    vector_list.pop_back();
    printer.PushAttribute( "Vectors", vector_list.c_str() );
  }
  if (scalar_list.size()>1) {
    scalar_list.pop_back();
    printer.PushAttribute( "Scalars", scalar_list.c_str() );
  }

  if (has_strengths) {
    printer.OpenElement( "DataArray" );
    printer.PushAttribute( "Name", "circulation" );
    printer.PushAttribute( "type", "Float32" );
    write_DataArray (printer, pts.get_str(), compress, asbase64);
    printer.CloseElement();	// DataArray
  }

  if (has_radii) {
    printer.OpenElement( "DataArray" );
    printer.PushAttribute( "Name", "radius" );
    printer.PushAttribute( "type", "Float32" );
    write_DataArray (printer, pts.get_rad(), compress, asbase64);
    printer.CloseElement();	// DataArray
  }

  printer.OpenElement( "DataArray" );
  printer.PushAttribute( "NumberOfComponents", "3" );
  printer.PushAttribute( "Name", "velocity" );
  printer.PushAttribute( "type", "Float32" );
  write_DataArray (printer, pts.get_vel(), compress, asbase64);
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

  std::cout << "Wrote " << pts.get_n() << " points to " << vtkfn.str() << std::endl;
  return vtkfn.str();
}


//
// write surface/panel data to a .vtu file
//
template <class S>
std::string write_vtu_panels(Surfaces<S> const& surf, const size_t file_idx,
                             const size_t frameno, const double time) {

  assert(surf.get_npanels() > 0 && "Inside write_vtu_panels with no panels");

  const bool compress = false;
  const bool asbase64 = true;

  bool has_vort_str = false;
  bool has_src_str = false;
  std::string prefix = "panel_";
  if (not surf.is_inert()) {
    has_vort_str = true;
    has_src_str = surf.have_src_str();
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

  // must choose one of these two formats
  printer.OpenElement( "UnstructuredGrid" );
  //printer.OpenElement( "PolyData" );

  // include simulation time here
  printer.OpenElement( "FieldData" );
  printer.OpenElement( "DataArray" );
  printer.PushAttribute( "type", "Float64" );
  printer.PushAttribute( "Name", "TimeValue" );
  printer.PushAttribute( "NumberOfTuples", "1" );
  {
    Vector<double> time_vec = {time};
    write_DataArray (printer, time_vec, false, false);
  }
  printer.CloseElement();	// DataArray
  printer.CloseElement();	// FieldData

  printer.OpenElement( "Piece" );
  printer.PushAttribute( "NumberOfPoints", std::to_string(surf.get_n()).c_str() );
  printer.PushAttribute( "NumberOfCells", std::to_string(surf.get_npanels()).c_str() );

  printer.OpenElement( "Points" );
  printer.OpenElement( "DataArray" );
  printer.PushAttribute( "NumberOfComponents", "3" );	// SEE THIS 3?!?!? It needs to be there.
  printer.PushAttribute( "Name", "position" );
  printer.PushAttribute( "type", "Float32" );
  write_DataArray (printer, surf.get_pos(), compress, asbase64);
  printer.CloseElement();	// DataArray
  printer.CloseElement();	// Points

  printer.OpenElement( "Cells" );

  printer.OpenElement( "DataArray" );
  printer.PushAttribute( "Name", "connectivity" );
  if (surf.get_n() <= std::numeric_limits<uint16_t>::max()) {
    printer.PushAttribute( "type", "UInt16" );
    std::vector<Int> const & idx = surf.get_idx();
    Vector<uint16_t> v(std::begin(idx), std::end(idx));
    write_DataArray (printer, v, compress, asbase64);
  } else {
    printer.PushAttribute( "type", "UInt32" );
    std::vector<Int> const & idx = surf.get_idx();
    Vector<uint32_t> v(std::begin(idx), std::end(idx));
    write_DataArray (printer, v, compress, asbase64);
  }
  printer.CloseElement();	// DataArray

  printer.OpenElement( "DataArray" );
  printer.PushAttribute( "Name", "offsets" );
  if (surf.get_n() <= std::numeric_limits<uint16_t>::max()) {
    printer.PushAttribute( "type", "UInt16" );
    Vector<uint16_t> v(surf.get_npanels());
    std::iota(v.begin(), v.end(), 1);
    std::transform(v.begin(), v.end(), v.begin(),
                   std::bind(std::multiplies<uint16_t>(), std::placeholders::_1, 2));
    write_DataArray (printer, v, compress, asbase64);
  } else {
    printer.PushAttribute( "type", "UInt32" );
    Vector<uint32_t> v(surf.get_npanels());
    std::iota(v.begin(), v.end(), 1);
    std::transform(v.begin(), v.end(), v.begin(),
                   std::bind(std::multiplies<uint32_t>(), std::placeholders::_1, 2));
    write_DataArray (printer, v, compress, asbase64);
  }
  printer.CloseElement();	// DataArray

  printer.OpenElement( "DataArray" );
  printer.PushAttribute( "Name", "types" );
  printer.PushAttribute( "type", "UInt8" );
  Vector<uint8_t> v(surf.get_npanels());
  std::fill(v.begin(), v.end(), 3);
  write_DataArray (printer, v, compress, asbase64);
  printer.CloseElement();	// DataArray

  printer.CloseElement();	// Cells


  printer.OpenElement( "CellData" );
  std::string vector_list;
  std::string scalar_list;

  if (has_vort_str) scalar_list.append("vortex sheet strength,");
  if (has_src_str) scalar_list.append("source sheet strength,");
  //if (has_radii) scalar_list.append("area,");
  vector_list.append("velocity,");

  if (vector_list.size()>1) {
    vector_list.pop_back();
    printer.PushAttribute( "Vectors", vector_list.c_str() );
  }
  if (scalar_list.size()>1) {
    scalar_list.pop_back();
    printer.PushAttribute( "Scalars", scalar_list.c_str() );
  }

  if (has_vort_str) {
    // put sheet strengths in here
    printer.OpenElement( "DataArray" );
    printer.PushAttribute( "Name", "vortex sheet strength" );
    printer.PushAttribute( "type", "Float32" );
    write_DataArray (printer, surf.get_vort_str(), compress, asbase64);
    printer.CloseElement();	// DataArray
  }

  if (has_src_str) {
    // put sheet strengths in here
    printer.OpenElement( "DataArray" );
    printer.PushAttribute( "Name", "source sheet strength" );
    printer.PushAttribute( "type", "Float32" );
    write_DataArray (printer, surf.get_src_str(), compress, asbase64);
    printer.CloseElement();	// DataArray
  }

  printer.OpenElement( "DataArray" );
  printer.PushAttribute( "NumberOfComponents", "3" );
  printer.PushAttribute( "Name", "velocity" );
  printer.PushAttribute( "type", "Float32" );
  write_DataArray (printer, surf.get_vel(), compress, asbase64);
  printer.CloseElement();	// DataArray

  printer.CloseElement();	// CellData


  printer.CloseElement();	// Piece
  printer.CloseElement();	// PolyData or UnstructuredGrid
  printer.CloseElement();	// VTKFile

  std::fclose(fp);

  std::cout << "Wrote " << surf.get_npanels() << " panels to " << vtkfn.str() << std::endl;
  return vtkfn.str();
}


//
// write a collection
//
template <class S>
void write_vtk_files(std::vector<Collection> const& coll, const size_t _index, const double _time,
                     std::vector<std::string>& _files) {

  size_t idx = 0;
  for (auto &elem : coll) {
    // is there a way to simplify this code with a lambda?
    //std::visit([=](auto& elem) { elem.write_vtk(); }, coll);

    // split on collection type
    if (std::holds_alternative<Points<S>>(elem)) {
      Points<S> const & pts = std::get<Points<S>>(elem);
      if (pts.get_n() > 0) {
        _files.emplace_back(write_vtu_points<S>(pts, idx++, _index, _time));
      }
    } else if (std::holds_alternative<Surfaces<S>>(elem)) {
      Surfaces<S> const & surf = std::get<Surfaces<S>>(elem);
      if (surf.get_npanels() > 0) {
        _files.emplace_back(write_vtu_panels<S>(surf, idx++, _index, _time));
      }
    }
  }
}

