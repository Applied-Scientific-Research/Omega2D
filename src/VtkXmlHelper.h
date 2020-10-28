/*
 * VtkXmlHelper.h - Write an XML-format VTK data file using TinyXML2
 *
 * (c)2019-20 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
 *            Blake B Hillier <blakehillier@mac.com>
 */

#pragma once

#include "Omega2D.h"
#include "VectorHelper.h"
#include "Collection.h"
#include "Points.h"
#include "Surfaces.h"
#include "Volumes.h"
#include "tinyxml2.h"
#include "VtkXmlWriter.h"
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

// write point data to a .vtu file
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
  VtkXmlWriter ptsWriter = VtkXmlWriter(vtkfn.str(), asbase64);
  // push comment with sim time?

  // include simulation time here
  ptsWriter.addElement("FieldData");

  {
    std::map<std::string, std::string> attribs = {{"type",           "Float64"},
                                                  {"Name",           "TimeValue"},
                                                  {"NumberOfTuples", "1"}};
    ptsWriter.addElement("DataArray", attribs);
    Vector<double> time_vec = {time};
    ptsWriter.writeDataArray(time_vec);
    // DataArray
    ptsWriter.closeElement();
  }
  // FieldData
  ptsWriter.closeElement();

  {
    std::map<std::string, std::string> attribs = {{"NumberOfPoints", std::to_string(pts.get_n()).c_str()},
                                                  {"NumberOfCells", std::to_string(pts.get_n()).c_str()}};
    ptsWriter.addElement("Piece", attribs);
  }
  
  ptsWriter.addElement("Points");
  {
    std::map<std::string, std::string> attribs = {{"NumberOfComponents", "3"},
                                                  {"Name",               "position"},
                                                  {"type",               "Float32"}};
    ptsWriter.addElement("DataArray", attribs);
    Vector<float> pos = ptsWriter.unpackArray(pts.get_pos());
    ptsWriter.writeDataArray(pos);
    // DataArray
    ptsWriter.closeElement();
  }
  // Points
  ptsWriter.closeElement();

  ptsWriter.addElement("Cells");
  
  // https://discourse.paraview.org/t/cannot-open-vtu-files-with-paraview-5-8/3759
  // apparently the Vtk format documents indicate that connectivities and offsets
  //   must be in Int32, not UIntAnything. Okay...
  {
    std::map<std::string, std::string> attribs = {{"Name", "connectivity"},
                                                  {"type", "Int32"}};
    ptsWriter.addElement("DataArray", attribs);
    Vector<int32_t> v(pts.get_n());
    std::iota(v.begin(), v.end(), 0);
    ptsWriter.writeDataArray(v);
    // DataArray
    ptsWriter.closeElement();
  }

  {
    std::map<std::string, std::string> attribs = {{"Name", "offsets"},
                                                  {"type", "Int32"}};
    ptsWriter.addElement("DataArray", attribs);
    Vector<int32_t> v(pts.get_n());
    std::iota(v.begin(), v.end(), 1);
    ptsWriter.writeDataArray(v);
    // DataArray
    ptsWriter.closeElement();
  }

  // except these, they can be chars
  {
    std::map<std::string, std::string> attribs = {{"Name", "types"},
                                                  {"type", "UInt8"}};
    ptsWriter.addElement("DataArray", attribs);
    Vector<uint8_t> v(pts.get_n());
    std::fill(v.begin(), v.end(), 1);
    ptsWriter.writeDataArray(v);
    // DataArray
    ptsWriter.closeElement();
  }
  // Cells
  ptsWriter.closeElement();

  {
    std::map<std::string, std::string> attribs = {{"Vectors", "velocity"}};
    std::string scalar_list;
    if (has_strengths) scalar_list.append("circulation,");
    if (has_radii) scalar_list.append("radius,");
    if (scalar_list.size()>1) {
      scalar_list.pop_back();
      attribs.insert({"Scalars", scalar_list});
    }
    ptsWriter.addElement("PointData", attribs);
  }

  if (has_strengths) {
    std::map<std::string, std::string> attribs = {{"Name", "circulation"},
                                                  {"type", "Float32"}};
    ptsWriter.addElement("DataArray", attribs);
    ptsWriter.writeDataArray(pts.get_str());
    // DataArray
    ptsWriter.closeElement();
  }

  if (has_radii) {
    std::map<std::string, std::string> attribs = {{"Name", "radius"},
                                                  {"type", "Float32"}};
    ptsWriter.addElement("DataArray", attribs);
    ptsWriter.writeDataArray(pts.get_rad());
    // DataArray
    ptsWriter.closeElement();
  }

  {
    std::map<std::string, std::string> attribs = {{"NumberOfComponents", "3"},
                                                  {"Name",               "velocity"},
                                                  {"type",               "Float32"}};
    ptsWriter.addElement("DataArray", attribs);
    Vector<float> vel = ptsWriter.unpackArray(pts.get_vel());
    ptsWriter.writeDataArray(vel);
  }
  // DataArray
  ptsWriter.closeElement();

  // Point Data 
  ptsWriter.closeElement();

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

  // Piece 
  ptsWriter.closeElement();
  ptsWriter.finish();
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
  VtkXmlWriter panelWriter = VtkXmlWriter(vtkfn.str(), asbase64);
  // push comment with sim time?

  // include simulation time here
  panelWriter.addElement("FieldData");
  {
    std::map<std::string, std::string> attribs = {{"type",           "Float64"},
                                                  {"Name",           "TimeValue"},
                                                  {"NumberOfTuples", "1"}};
    panelWriter.addElement("DataArray", attribs);
    Vector<double> time_vec = {time};
    panelWriter.writeDataArray(time_vec);
    panelWriter.closeElement();
  }
  // FieldData
  panelWriter.closeElement();

  {
    std::map<std::string, std::string> attribs = {{"NumberOfPoints", std::to_string(surf.get_n()).c_str()},
                                                  {"NumberOfCells", std::to_string(surf.get_npanels()).c_str()}};
    panelWriter.addElement("Piece", attribs);
  }

  panelWriter.addElement("Points");
  {
    std::map<std::string, std::string> attribs = {{"NumberOfComponents", "3"},
                                                  {"Name",               "position"},
                                                  {"type",               "Float32"}}; 
    panelWriter.addElement("DataArray", attribs);
    Vector<float> pos = panelWriter.unpackArray(surf.get_pos());
    panelWriter.writeDataArray(pos);
    panelWriter.closeElement();
  }
  // Points
  panelWriter.closeElement();

  panelWriter.addElement("Cells");
  // again, all connectivities and offsets must be Int32!
  {
    std::map<std::string, std::string> attribs = {{"Name", "connectivity"},
                                                  {"type", "Int32"}};
    panelWriter.addElement("DataArray", attribs);
    std::vector<Int> const & idx = surf.get_idx();
    Vector<int32_t> v(std::begin(idx), std::end(idx));
    panelWriter.writeDataArray(v);
    panelWriter.closeElement();
  }

  {
    std::map<std::string, std::string> attribs = {{"Name", "offsets"},
                                                  {"type", "Int32"}};
    panelWriter.addElement("DataArray", attribs);
    Vector<int32_t> v(surf.get_npanels());
    std::iota(v.begin(), v.end(), 1);
    std::transform(v.begin(), v.end(), v.begin(),
                   std::bind(std::multiplies<int32_t>(), std::placeholders::_1, 2));
    panelWriter.writeDataArray(v);
    panelWriter.closeElement();
  }

  {
    std::map<std::string, std::string> attribs = {{"Name", "types"},
                                                  {"type", "UInt8"}};
    panelWriter.addElement("DataArray", attribs);
    Vector<uint8_t> v(surf.get_npanels());
    std::fill(v.begin(), v.end(), 3);
    panelWriter.writeDataArray(v);
    panelWriter.closeElement();
  }
  // Cells
  panelWriter.closeElement();

  {
    std::map<std::string, std::string> attribs = {{"Vectors", "velocity"}};
    std::string scalar_list;
    if (has_vort_str) scalar_list.append("vortex sheet strength,");
    if (has_src_str) scalar_list.append("source sheet strength,");
    //if (has_radii) scalar_list.append("area,");
    if (scalar_list.size()>1) {
      scalar_list.pop_back();
      attribs.insert({"Scalars", scalar_list});
    }
    panelWriter.addElement("CellData", attribs);
  }

  if (has_vort_str) {
    std::map<std::string, std::string> attribs = {{"Name", "vortex sheet strength"},
                                                  {"type", "Float32"}};
    panelWriter.addElement("DataArray", attribs);
    panelWriter.writeDataArray(surf.get_vort_str());
    panelWriter.closeElement();
  }

  if (has_src_str) {
    std::map<std::string, std::string> attribs = {{"Name", "source sheet strength"},
                                                  {"type", "Float32"}};
    panelWriter.addElement("DataArray", attribs);
    panelWriter.writeDataArray(surf.get_src_str());
    panelWriter.closeElement();
  }

  {
    std::map<std::string, std::string> attribs = {{"NumberOfComponents", "3"},
                                                  {"Name",               "velocity"},
                                                  {"type",               "Float32"}}; 
    panelWriter.addElement("DataArray", attribs);
    Vector<float> vel = panelWriter.unpackArray(surf.get_vel());
    panelWriter.writeDataArray(vel);
    panelWriter.closeElement();
  }
  // CellData
  panelWriter.closeElement();
  // Piece 
  panelWriter.closeElement();

  panelWriter.finish();
  std::cout << "Wrote " << surf.get_npanels() << " panels to " << vtkfn.str() << std::endl;
  return vtkfn.str();
}

// write grid data to a .vtk file
template <class S>
std::string write_vtk_grid(Volumes<S> const& grid, const size_t file_idx,
                           const size_t frameno, const double time) {

  assert(grid.get_nelems() > 0 && "Inside write_vtk_grid with no elements");
  const bool asbase64 = true;

  // generate file name
  std::string prefix = "grid_";
  std::stringstream vtkfn;
  vtkfn << prefix << std::setfill('0') << std::setw(2) << file_idx << "_" << std::setw(5) << frameno << ".vtu";
  VtkXmlWriter gridWriter = VtkXmlWriter(vtkfn.str(), asbase64);

  // include simulation time here
  gridWriter.addElement("FieldData");
  {
    std::map<std::string, std::string> attribs = {{"type", "Float64"},
                                                  {"Name", "TimeValue"},
                                                  {"NumberOfTuples", "1"}};
    gridWriter.addElement("DataArray", attribs);
    Vector<double> time_vec = {time};
    gridWriter.writeDataArray(time_vec);
    gridWriter.closeElement();
  }
  // FieldData
  gridWriter.closeElement();

  {
    std::map<std::string, std::string> attribs = {{"NumberOfPoints", std::to_string(grid.get_n())},
                                                  {"NumberOfCells", std::to_string(grid.get_nelems())}};
    gridWriter.addElement("Piece", attribs);
  }

  gridWriter.addElement("Points");
  {
    std::map<std::string, std::string> attribs = {{"NumberOfComponents", "3"},
                                                  {"Name", "position"},
                                                  {"type", "Float32"}};
    gridWriter.addElement("DataArray", attribs);
    Vector<float> pos = gridWriter.unpackArray(grid.get_pos());
    gridWriter.writeDataArray(pos);
    gridWriter.closeElement();
  }
  // Points
  gridWriter.closeElement();

  gridWriter.addElement("Cells");

  // again, all connectivities and offsets must be Int32!
  {
    std::map<std::string, std::string> attribs = {{"Name", "connectivity"},
                                                  {"type", "Int32"}};
    gridWriter.addElement("DataArray", attribs);
    std::vector<Int> const & idx = grid.get_idx();
    Vector<int32_t> v(std::begin(idx), std::end(idx));
    gridWriter.writeDataArray(v);
    // DataArray
    gridWriter.closeElement();
  }

  {
    std::map<std::string, std::string> attribs = {{"Name", "offsets"},
                                                  {"type", "Int32"}};
    gridWriter.addElement("DataArray", attribs);
    Vector<int32_t> v(grid.get_nelems());
    // vector of 1 to n
    std::iota(v.begin(), v.end(), 1);
    std::transform(v.begin(), v.end(), v.begin(),
                   std::bind(std::multiplies<int32_t>(), std::placeholders::_1, 4));
    gridWriter.writeDataArray(v);
    // DataArray
    gridWriter.closeElement();
  }

  {
    std::map<std::string, std::string> attribs = {{"Name", "types"},
                                                  {"type", "UInt8"}};
    gridWriter.addElement("DataArray", attribs);
    Vector<uint8_t> v(grid.get_nelems());
    std::fill(v.begin(), v.end(), 9);
    gridWriter.writeDataArray(v);
    // DataArray
    gridWriter.closeElement();
  }
  // Cells
  gridWriter.closeElement();
  gridWriter.addElement("PointData");
  
  {
    std::map<std::string, std::string> attribs = {{"NumberOfComponents", "3"},
                                                  {"Name", "velocity"},
                                                  {"type", "Float32"}};
    gridWriter.addElement("DataArray", attribs);
    Vector<float> vel = gridWriter.unpackArray(grid.get_vel());
    gridWriter.writeDataArray(vel);
    // DataArray
    gridWriter.closeElement();
  }
 
  // PointData
  gridWriter.closeElement();
  // Piece
  gridWriter.closeElement();

  gridWriter.finish();
  std::cout << "Wrote " << grid.get_nelems() << " elements to " << vtkfn.str() << std::endl;
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
    } else if (std::holds_alternative<Volumes<S>>(elem)) {
      Volumes<S> const & cells = std::get<Volumes<S>>(elem);
      if (cells.get_nelems() > 0) {
        _files.emplace_back(write_vtk_grid<S>(cells, idx++, _index, _time));
      }
    }
  }
}

