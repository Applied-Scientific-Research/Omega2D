/*
 * VtkXmlWriter.h - Write an XML-format VTK data file using TinyXML2
 *
 * This class manages a writer object for vtk using the xml format.
 * Currently only supports 64 encoding and unstructured grids.
 * The next step is to start generalizing bigger writing processes for the vtu format.
 * 
 * (c)2019-20 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
 *            Blake B Hillier <blakehillier@mac.com>
 */

#pragma once

#include "cppcodec/base64_rfc4648.hpp"
#include "tinyxml2/tinyxml2.h"

#include <cstdint>
#include <cstdio>
#include <map>
#include <string>
#include <vector>


class VtkXmlWriter {
public:
  VtkXmlWriter(const std::string _file, const bool _base64) {
    // prepare file pointer and p
    m_fp = std::fopen(_file.c_str(), "wb");
    m_p = std::make_unique<tinyxml2::XMLPrinter>(m_fp);

    // Write document descriptors 
    m_p->PushHeader(false, true);
    m_p->OpenElement( "VTKFile" );
    m_p->PushAttribute( "type", "UnstructuredGrid" );
    m_p->PushAttribute( "version", "0.1" );
    m_p->PushAttribute( "byte_order", "LittleEndian" );
    // note this is still unsigned even though all indices later are signed!
    m_p->PushAttribute( "header_type", "UInt32" );
    m_p->OpenElement( "UnstructuredGrid" );

    // Set writer properties 
    m_base64 = _base64;
    m_openElements = 2;
  }

  void addElement(const char* _text, std::map<std::string, std::string> params=std::map<std::string, std::string>()) {
    m_p->OpenElement(_text);
    if (params.size()>0) {
      for (auto it=params.begin(); it !=params.end(); ++it) {
        m_p->PushAttribute(it->first.c_str(), it->second.c_str());
      }
    }
    m_openElements++;
  }

  // interleave arrays and write to the vtk file
  template <class S, long unsigned int I>
  Vector<float> unpackArray(std::array<Vector<S>,I> const & _data) {
    // interleave the (2,3) vectors into a new one
    { assert((I==2 || I==3) && "ERROR with template I"); }
    Vector<float> newvec;
    newvec.resize(3 * _data[0].size());
    for (size_t i=0; i<_data[0].size(); ++i) {
      newvec[3*i+0] = _data[0][i];
      newvec[3*i+1] = _data[1][i];
      if (I == 2) {
        newvec[3*i+2] = 0.0;
      } else if (I == 3) {
        newvec[3*i+2] = _data[2][i];
      }
    }
    return newvec;
  }

  template <class S>
  void writeDataArray (Vector<S> const & _data) {
    using base64 = cppcodec::base64_rfc4648;
  
    // why would you ever want to use base64 for floats and such? so wasteful.
    if (m_base64) {
      m_p->PushAttribute( "format", "binary" );
      std::string encoded = base64::encode((const char*)_data.data(), (size_t)(sizeof(_data[0])*_data.size()));
  
      // encode and write the UInt32 length indicator - the length of the base64-encoded blob
      uint32_t encoded_length = encoded.size();
      std::string header = base64::encode((const char*)(&encoded_length), (size_t)(sizeof(uint32_t)));
      //std::cout << "base64 length of header: " << header.size() << std::endl;
      //std::cout << "base64 length of data: " << encoded.size() << std::endl;
  
      m_p->PushText( " " );
      m_p->PushText( header.c_str() );
      m_p->PushText( encoded.c_str() );
      m_p->PushText( " " );
  
      if (false) {
        std::vector<uint8_t> decoded = base64::decode(encoded);
        S* decoded_floats = reinterpret_cast<S*>(decoded.data());
    
        std::cout << "Encoding an array of " << _data.size() << " floats to: " << encoded << std::endl;
        std::cout << "Decoding back into array of floats, starting with: " << decoded_floats[0] << " " << decoded_floats[1] << std::endl;
        //std::cout << "Decoding back into array of " << decoded.size() << " floats, starting with: " << decoded[0] << std::endl;
      }
    } else {
      m_p->PushAttribute( "format", "ascii" );
  
      // write the data
      m_p->PushText( " " );
      for (size_t i=0; i<_data.size(); ++i) {
        m_p->PushText( _data[i] );
        m_p->PushText( " " );
      }
    }
  }

  void closeElement() {
    m_p->CloseElement();
    m_openElements--;
  }

  void finish() {
    if (m_openElements != 2) { std::cout << "WARNING: " << m_openElements << " WERE NOT CLOSED" << std::endl; }
    while (m_openElements) {
      this->closeElement();
    }
    std::fclose(m_fp);
  }

private:
  bool m_base64;
  std::FILE* m_fp;
  int m_openElements;
  std::unique_ptr<tinyxml2::XMLPrinter> m_p;
};

