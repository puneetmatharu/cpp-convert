// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2022 Matthias Heil and Andrew Hazel
// LIC//
// LIC// This library is free software; you can redistribute it and/or
// LIC// modify it under the terms of the GNU Lesser General Public
// LIC// License as published by the Free Software Foundation; either
// LIC// version 2.1 of the License, or (at your option) any later version.
// LIC//
// LIC// This library is distributed in the hope that it will be useful,
// LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
// LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// LIC// Lesser General Public License for more details.
// LIC//
// LIC// You should have received a copy of the GNU Lesser General Public
// LIC// License along with this library; if not, write to the Free Software
// LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// LIC// 02110-1301  USA.
// LIC//
// LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
// LIC//
// LIC//====================================================================
//
// Copyright (C) 2007- Angelo Simone (TU Delft)
// Licensed under the GNU LGPL Version 2.1
//
// Script for converting from oomph-lib tecplot format to VTK XML
//
// Report bugs to a.simone@tudelft.nl
//
// The script can convert
// - 2d and 3d quad and tetra meshes
// - 2d and 3d points from all meshes
//
// The script can NOT convert
// - 1d meshes
//
// Scan the script for "Assumption" and "Note"
//==========================================================================

// TODO: Test loading with Python Paraview!
// TODO: Add functionality to save in binary format
// TODO: Run iwyu

#include <array>
#include <chrono>
#include <cstdint>
#include <stdexcept>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <optional>
#include <regex>
#include <vector>

#include "oomph_utilities.h"

namespace fs = std::filesystem;

/************************************************************************************
 * @brief
 *
 ************************************************************************************/
void print_usage()
{
  std::cout << R"(

NAME
        oomph-convert.py - script for converting from oomph-lib tecplot format to VTK XML.


SYNOPSIS
        oomph-convert.py [OPTION] [-i input_file.dat] [-o output_file]


OPTIONS
        -p2 or -p3 outputs only points in 2D or 3D (.vtp)
        -h         display this help text and exit
        -z         add trailing zeros to the output filename
        -f         force overwrite existing files

        By default output files will overwrite old ones only if the input file is
        newer than the output file.


TYPICAL USAGE EXAMPLES
        oomph-convert.py -h                             -> display help text
        oomph-convert.py soln12.dat soln12.vtu          -> generate soln12.vtu
        oomph-convert.py -z soln12.dat soln12.vtu       -> generate soln00012.vtu
        oomph-convert.py soln12.dat                     -> generate soln12.vtu
        oomph-convert.py -z soln12.dat                  -> generate soln00012.vtu
        oomph-convert.py soln12.dat nsol2.vtu           -> generate nsol2.vtu
        oomph-convert.py -z soln12.dat nsol2.vtu        -> generate nsol00002.vtu
        oomph-convert.py soln.dat                       -> generate soln.vtu
        oomph-convert.py -z soln.dat                    -> generate soln00000.vtu
        oomph-convert.py -p3 soln12.dat soln12.vtp          -> generate soln12.vtp
        oomph-convert.py -p2 -z soln12.dat soln12.vtp       -> generate soln00012.vtp
        oomph-convert.py -p3 soln12.dat                     -> generate soln12.vtp
        oomph-convert.py -p2 -z soln12.dat                  -> generate soln00012.vtp
        oomph-convert.py -p3 soln12.dat nsol2.vtp           -> generate nsol2.vtp
        oomph-convert.py -p2 -z soln12.dat nsol2.vtp        -> generate nsol00002.vtp
        oomph-convert.py -p3 soln.dat                       -> generate soln.vtp
        oomph-convert.py -p2 -z soln.dat                    -> generate soln00000.vtp
    )" << std::endl;
}


/************************************************************************************
 * @brief
 *
 ************************************************************************************/
namespace StringHelpers
{
  const std::string Whitespace = " \t\n\r\f\v";

  /************************************************************************************
   * @brief trim from end of string (right)
   *
   * @param s:
   * @param t:
   * @return std::string&:
   ************************************************************************************/
  inline std::string rtrim(std::string s, std::string t = Whitespace)
  {
    s.erase(s.find_last_not_of(t) + 1);
    return s;
  }

  /************************************************************************************
   * @brief trim from beginning of string (left)
   *
   * @param s:
   * @param t:
   * @return std::string&:
   ************************************************************************************/
  inline std::string ltrim(std::string s, std::string t = Whitespace)
  {
    s.erase(0, s.find_first_not_of(t));
    return s;
  }

  /************************************************************************************
   * @brief trim from both ends of string (right then left)
   *
   * @param s:
   * @param t:
   * @return std::string&:
   ************************************************************************************/
  inline std::string trim(std::string s, std::string t = Whitespace)
  {
    return ltrim(rtrim(s, t), t);
  }

  /************************************************************************************
   * @brief
   *
   * @param str
   ************************************************************************************/
  std::vector<std::string> split_string(const std::string& text,
                                        const char& delimeter = ' ')
  {
    std::stringstream text_stream(text);
    std::string word;
    std::vector<std::string> split_string;
    while (std::getline(text_stream, word, delimeter))
    {
      split_string.push_back(word);
    }
    return split_string;
  }

  /************************************************************************************
   * @brief
   *
   ************************************************************************************/
  std::string repeat(const std::string& s, unsigned n_repeat)
  {
    std::string output = "";
    for (unsigned i = 0; i < n_repeat; i++) output += s;
    return output;
  }

  /************************************************************************************
   * @brief
   *
   * @param text
   * @return std::string
   ************************************************************************************/
  std::string lower(const std::string& text)
  {
    std::string text_copy = text;
    std::transform(text_copy.begin(),
                   text_copy.end(),
                   text_copy.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return text_copy;
  }

  /************************************************************************************
   * @brief
   *
   ************************************************************************************/
  double string_to_double(const std::string& text)
  {
    return std::stod(std::regex_replace(lower(text), std::regex("d"), "e"));
  };
} // namespace StringHelpers


/************************************************************************************
 * @brief This class acts as a namespace for the definition Vtk XML headers and
 * footers.
 ************************************************************************************/
namespace VtpXml
{
  // clang-format off
  const std::string Header{"<?xml version=\"1.0\"?>\n<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">"};
  const std::string PolyDataHeader{R"(<PolyData>)"};
  auto PieceHeader = [](const uint64_t& number_of_points,
                        const uint64_t& number_of_lines) -> const std::string
  {
    return "<Piece NumberOfPoints=\"" + std::to_string(number_of_points) +
           "\" NumberOfVerts=\"0\" NumberOfLines=\"" +
           std::to_string(number_of_lines) +
           "\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">";
  };
  const std::string PointsHeader{"<Points>\n<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">"};
  const std::string PointsFooter{"</DataArray>\n</Points>"};
  const std::string PointDataHeader{R"(<PointData>)"};
  auto FieldHeader = [](const unsigned& field) -> const std::string
  {
    return "<DataArray type=\"Float32\" Name=\"V" + std::to_string(field) +
           "\" format=\"ascii\">";
  };
  const std::string FieldFooter{R"(</DataArray>)"};
  const std::string PointDataFooter{R"(</PointData>)"};
  const std::string VertsHeader{R"(<Verts>)"};
  const std::string VertsFooter{R"(</Verts>)"};
  const std::string LinesHeader{R"(<Lines>)"};
  const std::string LinesFooter{R"(</Lines>)"};
  const std::string StripsHeader{R"(<Strips>)"};
  const std::string StripsFooter{R"(</Strips>)"};
  const std::string PolysHeader{R"(<Polys>)"};
  const std::string PolysFooter{R"(</Polys>)"};
  const std::string PieceFooter{R"(</Piece>)"};
  const std::string PolyDataFooter{R"(</PolyData>)"};
  const std::string Footer{R"(</VTKFile>)"};
  // PM: Unused
  // const std::string CellDataHeader{R"(<CellData>)"};
  // const std::string CellDataFooter{R"(</CellData>)"};
  // clang-format on
}; // namespace VtpXml


/************************************************************************************
 * @brief This class acts as a namespace for the definition Vtk XML headers and
 * footers.
 ************************************************************************************/
namespace VtkXml
{
  // clang-format off
  const std::string Header{"<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">"};
  const std::string UnstructuredGridHeader{R"(<UnstructuredGrid>)"};
  auto PieceHeader = [](const uint64_t& number_of_points,
                        const uint64_t& number_of_cells) -> const std::string
  {
    return "<Piece NumberOfPoints=\"" + std::to_string(number_of_points) +
           "\" NumberOfCells=\"" + std::to_string(number_of_cells) + "\">";
  };
  const std::string PointsHeader{"<Points>\n<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">"};
  const std::string PointsFooter{"</DataArray>\n</Points>"};
  const std::string CellsHeader{R"(<Cells>)"};
  const std::string ConnectivityHeader{R"(<DataArray type="Int32" Name="connectivity" format="ascii">)"};
  const std::string ConnectivityFooter{R"(</DataArray>)"};
  const std::string OffsetsHeader{R"(<DataArray type="Int32" Name="offsets" format="ascii">)"};
  const std::string OffsetsFooter{R"(</DataArray>)"};
  const std::string TypesHeader{R"(<DataArray type="UInt8" Name="types" format="ascii">)"};
  const std::string TypesFooter{R"(</DataArray>)"};
  const std::string CellsFooter{R"(</Cells>)"};
  const std::string PointDataHeader{R"(<PointData>)"};
  auto FieldHeader = [](const unsigned& field) -> const std::string
  {
    return "<DataArray type=\"Float32\" Name=\"V" + std::to_string(field) +
           "\" format=\"ascii\">";
  };
  const std::string FieldFooter{R"(</DataArray>)"};
  const std::string PointDataFooter{R"(</PointData>)"};
  const std::string PieceFooter{R"(</Piece>)"};
  const std::string UnstructuredGridFooter{R"(</UnstructuredGrid>)"};
  const std::string Footer{R"(</VTKFile>)"};
  // clang-format on
}; // namespace VtkXml


namespace Tecplot
{
  /************************************************************************************
   * @brief This class handles Tecplot node informations.
   *
   ************************************************************************************/
  struct TecplotNode
  {
    // Node coordinates
    std::array<float, 3> coordinates{0.0, 0.0, 0.0};

    // Fields values
    std::vector<float> Fields{};
  };

  /************************************************************************************
   * @brief
   *
   ************************************************************************************/
  enum ZoneFormat
  {
    TriangleOrTetrahedron = 1,
    Other = 2,
  };

  /************************************************************************************
   * @brief
   *
   ************************************************************************************/
  enum VTKCellType
  {
    Vertex = 1,
    PolyVertex = 2,
    Line = 3,
    PolyLine = 4,
    Triangle = 5,
    TriangleStrip = 6,
    Polygon = 7,
    Pixel = 8,
    Quad = 9,
    Tetrahedron = 10,
    Voxel = 11,
    Hexahedron = 12,
    Wedge = 13,
    Pyramid = 14,
  };

  /************************************************************************************
   * @brief
   *
   ************************************************************************************/
  struct Connectivity
  {
    Connectivity(const std::vector<std::string>& connectivities)
    {
      Connectivity_list.resize(connectivities.size());
      std::transform(connectivities.begin(),
                     connectivities.end(),
                     Connectivity_list.begin(),
                     [](const std::string& s)
                     {
                       try
                       {
                         return std::stoi(s);
                       }
                       catch (std::invalid_argument)
                       {
                         std::cout << "\n\nUnable to convert argument '" << s
                                   << "' to an integer!\n"
                                   << std::endl;
                         exit(1);
                       }
                     });
    }

    /************************************************************************************
     * @brief
     *
     * @return std::vector<int32_t>::size_type:
     ************************************************************************************/
    std::vector<int32_t>::size_type size() const
    {
      return Connectivity_list.size();
    }

    /************************************************************************************
     * @brief
     *
     * @param idx:
     * @return double:
     ************************************************************************************/
    int32_t operator[](int idx) const
    {
#ifdef RANGE_CHECKING
      return Connectivity_list.at(idx);
#else
      return Connectivity_list[idx];
#endif
    }

    /************************************************************************************
     * @brief
     *
     * @param idx:
     * @return v_ref:
     ************************************************************************************/
    int32_t& operator[](int idx)
    {
#ifdef RANGE_CHECKING
      return Connectivity_list.at(idx);
#else
      return Connectivity_list[idx];
#endif
    }

    /************************************************************************************
     * @brief
     *
     * @return std::vector<int32_t>&:
     ************************************************************************************/
    const std::vector<int32_t>& data() const
    {
      return Connectivity_list;
    }

  private:
    std::vector<int32_t> Connectivity_list;
  };

  /************************************************************************************
   * @brief This class handles Tecplot zone informations.
   *
   ************************************************************************************/
  struct TecplotZone
  {
    std::vector<int32_t> Edges{};
    std::vector<TecplotNode> Nodes{};
    std::vector<uint8_t> Dimension{};
    std::vector<VTKCellType> Cell_type{};
    std::vector<ZoneFormat> Cell_format{};
    std::vector<uint32_t> N_cell{};
    std::vector<Connectivity> Connectivities{};

    /************************************************************************************
     * @brief
     *
     * @return int:
     ************************************************************************************/
    std::vector<TecplotNode>::size_type nnode() const
    {
      return Nodes.size();
    }
  };

  /************************************************************************************
   * @brief Parse the zone of the given file starting from current position.
   *
   * @param file_reader: A file open for reading.
   * @param line_number: The current position in the file as line number (used
   *        for error messages).
   * @return std::tuple<std::optional<TecplotZone>, uint64_t>: a TecplotZone
   *         object if a zone has been found before the end of the file,
   *         otherwise nullopt. The second output will be the update line
   *         position in the output file.
   ************************************************************************************/
  std::optional<TecplotZone> maybe_parse_next_tecplot_zone(
    std::ifstream& file_reader, uint64_t& line_number)
  {
    //-----------------------------------------------------------------------
    // Seek to the next Tecplot zone
    //-----------------------------------------------------------------------
    std::string header;
    while (true)
    {
      std::getline(file_reader, header);
      line_number += 1;

      // Have we reached the end of the file?
      if (!file_reader)
      {
        return std::nullopt;
      }

      // Check if this line describes a zone
      if ((header.length() > 3) && (header.substr(0, 4) == "ZONE"))
      {
        // We got a zone!
        break;
      }
    }

    //-----------------------------------------------------------------------
    // Create the zone
    //-----------------------------------------------------------------------
    TecplotZone zone{};

    // Read header informations. Strip the "ZONE" prefix and any surrounding
    // whitespaces, then split based on commas
    auto edges_info = header.substr(4);
    auto edges =
      StringHelpers::split_string(StringHelpers::trim(edges_info), ',');

    // Assumption: Only two zone formats are possible.
    // The first (format 1) is defined with a zone header of the type
    // ZONE N=15, E=16, F=FEPOINT, ET=TRIANGLE
    // In this case, edges = ['N=15', ' E=16', ' F=FEPOINT', '
    // ET=TRIANGLE'] This zone has been coded for two ET keywords:
    // TRIANGLE and TETRAHEDRON The second (format 2) is defined in terms
    // of the IJK indices and a typical zone reads ZONE I=5, J=5, K=5 One
    // index, I, indicates a line element, 2 indices, I and J, indicate a
    // 4 node quadrilateral element, and three indices, I J and K,
    // indicate a tetrahedral element.
    VTKCellType cell;
    uint8_t dim;
    ZoneFormat format;
    uint32_t n_node;
    uint32_t n_cell;

    if (header.find("TRIANGLE") != std::string::npos)
    {
      // Define data
      cell = VTKCellType::Triangle;
      dim = 2;
      format = ZoneFormat::TriangleOrTetrahedron;
      n_node =
        std::stoi(StringHelpers::trim(StringHelpers::trim(edges[0]), "N="));
      n_cell =
        std::stoi(StringHelpers::trim(StringHelpers::trim(edges[1]), "E="));

      // Append data to zone
      zone.Dimension.emplace_back(dim);
      zone.Cell_type.emplace_back(cell);
      zone.Cell_format.emplace_back(format);
      zone.N_cell.emplace_back(n_cell);
    }
    else if (header.find("TETRAHEDRON") != std::string::npos)
    {
      // Define data
      cell = VTKCellType::Tetrahedron; // VTK_TETRAHEDRON
      dim = 3;
      format = ZoneFormat::TriangleOrTetrahedron;
      n_node =
        std::stoi(StringHelpers::trim(StringHelpers::trim(edges[0]), "N="));
      n_cell =
        std::stoi(StringHelpers::trim(StringHelpers::trim(edges[1]), "E="));

      // Append data to zone
      zone.Dimension.emplace_back(dim);
      zone.Cell_type.emplace_back(cell);
      zone.Cell_format.emplace_back(format);
      zone.N_cell.emplace_back(n_cell);
    }
    else
    {
      size_t n_edge = edges.size();
      if ((n_edge < 1) || (n_edge > 3))
      {
        throw std::runtime_error("Parsing failed at line_number " +
                                 std::to_string(line_number) +
                                 ": [Invalid zone header] wrong edges count! "
                                 "Try to convert with -p option");
      }

      format = ZoneFormat::Other;

      // define dimension of zone cell by counting the IJK indices
      dim = n_edge;

      if (dim == 1) cell = VTKCellType::Line;
      if (dim == 2) cell = VTKCellType::Quad;
      if (dim == 3) cell = VTKCellType::Hexahedron;

      // append dimension and cell type/format to the zone
      zone.Dimension.emplace_back(dim);
      zone.Cell_type.emplace_back(cell);
      zone.Cell_format.emplace_back(format);

      // extract information from string
      std::array<char, 3> labels = {'I', 'J', 'K'};
      unsigned index{0};
      for (auto& edge : edges)
      {
        edge = StringHelpers::trim(edge);
        if ((edge.length() < 2) || (edge[0] != labels[index]))
        {
          throw std::runtime_error(
            "Parsing failed at line_number " + std::to_string(line_number) +
            ": [Invalid zone header] wrong edges format " + edges_info +
            ". Try to convert with -p option");
        }

        try
        {
          zone.Edges.emplace_back(std::stoi(edge.substr(2)));
        }
        catch (...)
        {
          throw std::runtime_error(
            "Parsing failed at line_number " + std::to_string(line_number) +
            ": [Invalid zone header] wrong edges format " + edges_info +
            ". Try to convert with -p option");
        }
        index++;
      }

      //-------------------------------------------------------------------
      // Node count in cell
      //-------------------------------------------------------------------
      n_node = 1;
      for (const auto& edge : zone.Edges)
      {
        n_node *= edge;
      }

      //-------------------------------------------------------------------
      // Cell count in zone
      //-------------------------------------------------------------------
      n_cell = 1;
      for (const auto& edge : zone.Edges)
      {
        n_cell *= (edge - 1);
      }
      zone.N_cell.emplace_back(n_cell);
    }

    //-----------------------------------------------------------------------
    // Parse nodes
    //-----------------------------------------------------------------------
    for (int n = 0; n < n_node; n++)
    {
      std::string data;
      std::getline(file_reader, data);
      line_number += 1;

      if (data.empty() || (StringHelpers::trim(data) == ""))
      {
        std::getline(file_reader, data); // ...in case of empty line
        line_number += 1;
      }

      auto data_segments =
        StringHelpers::split_string(StringHelpers::trim(data));

      if (data_segments.size() < dim)
      {
        throw std::runtime_error("Parsing failed at line_number " +
                                 std::to_string(line_number) +
                                 ": [Invalid zone] not enough values for "
                                 "this node! Try to convert with -p option");
      }

      TecplotNode node;
      unsigned i = 0;
      for (const auto& value : data_segments)
      {
        try
        {
          if (i < dim)
          {
            node.coordinates[i] = std::stod(value);
          }
          else
          {
            node.Fields.emplace_back(std::stod(value));
          }
        }
        catch (...)
        {
          throw std::runtime_error("Parsing failed at line_number " +
                                   std::to_string(line_number) +
                                   ": [Invalid zone] wrong node values! Try "
                                   "to convert with -p option");
        }
        i++;
      }

      // Append this node to the zone
      zone.Nodes.emplace_back(node);
    }

    //-----------------------------------------------------------------------
    // Parse connectivities (format == 1)
    //-----------------------------------------------------------------------
    if (format == ZoneFormat::TriangleOrTetrahedron)
    {
      for (int i = 0; i < n_cell; i++)
      {
        std::string data;
        std::getline(file_reader, data);
        line_number += 1;

        std::cout << "data: " << data << std::endl;
        if (data.empty() || (StringHelpers::trim(data) == ""))
        {
          std::getline(file_reader, data); // ...in case of empty line
          line_number += 1;
        }
        std::cout << "data: " << data << std::endl;

        auto data_segments =
          StringHelpers::split_string(StringHelpers::trim(data));

        if (cell == 3) n_node = 2; // VTK_LINE
        if (cell == 5) n_node = 3; // VTK_TRIANGLE
        if (cell == 9) n_node = 4; // VTK_QUAD
        if (cell == 10) n_node = 4; // VTK_TETRAHEDRON
        if (cell == 12) n_node = 8; // VTK_HEXAHEDRON

        if (data_segments.size() != n_node)
        {
          throw std::runtime_error(
            "Parsing failed at line_number " + std::to_string(line_number) +
            ": [Invalid zone] wrong connectivity list for this "
            "element!");
        }

        // Append these connectivities to the zone
        zone.Connectivities.emplace_back(data_segments);
      }
    }
    return zone;
  }

  /************************************************************************************
   * @brief
   *
   * @param file:
   * @param line:
   * @param dim:
   * @return std::tuple<TecplotNode, int>:
   ************************************************************************************/
  std::optional<TecplotNode> maybe_parse_next_tecplot_node(
    std::ifstream& file_reader, uint64_t& line_number, const unsigned& dim)
  {
    static std::regex float_point_pattern{
      R"([-+]?[0-9]*\.?[0-9]*([e,E][-+]?[0-9]+)?)"};

    //-----------------------------------------------------------------------
    // Create the point
    //-----------------------------------------------------------------------
    TecplotNode point{};

    //-----------------------------------------------------------------------
    // Seek to the next point
    //-----------------------------------------------------------------------
    while (true)
    {
      std::string line;
      std::getline(file_reader, line);
      line_number += 1;

      // We reach the end of the file
      if (!file_reader)
      {
        return std::nullopt;
      }

      // Boolean indicating that the line contains only numerical data
      bool is_fp_line{true};

      // Read the line and split into individual entries
      std::vector<std::string> values =
        StringHelpers::split_string(StringHelpers::trim(line));

      // Check that all entries are floating points numbers
      for (const auto& value : values)
      {
        // Define regular expression for floating point number
        std::smatch match_index;
        if (std::regex_match(value, match_index, float_point_pattern))
        {
          // If the matched string doesn't cover the entire number
          if (match_index.str(0).length() != value.length())
          {
            is_fp_line = false;
            break;
          }
        }
        else
        {
          is_fp_line = false;
          break;
        }
      }

      // If it's a line containing only floating point numbers:
      // Extract coordinates and values
      if (is_fp_line)
      {
        unsigned i = 0;
        for (const auto& value : values)
        {
          try
          {
            // Coordinate
            if (i < dim)
            {
              point.coordinates[i] = StringHelpers::string_to_double(value);
            }
            // Value
            else
            {
              point.Fields.emplace_back(StringHelpers::string_to_double(value));
            }
          }
          catch (...)
          {
            std::string error_message{"Parsing failed at line_number " +
                                      std::to_string(line_number) +
                                      ": [Invalid zone] wrong node values!"};
            throw std::invalid_argument(error_message);
          }

          i++;
        }
        return point;
      }
    }
  }

  /************************************************************************************
   * @brief Converts Tecplot file generates by oomph-lib into Vtk XML file.
   *
   * @param input_filename:
   * @param output_filename:
   ************************************************************************************/
  void tecplot_to_vtkxml(const fs::path& input_filename,
                         const fs::path& output_filename)
  {
    //---------------------------------------------------------------------------
    // Retrieve the points from the input Tecplot file
    //---------------------------------------------------------------------------
    std::ifstream reader(input_filename);
    if (!reader.is_open())
    {
      throw std::runtime_error("Failed to open input file for reading!");
    }

    std::cout << "Parse input file for Tecplot zones........" << std::flush;
    uint64_t line_num{0};
    int64_t n_ignored_line{0};
    std::vector<TecplotZone> zones{};

    while (true)
    {
      auto temp_line_num = line_num;
      auto zone = std::move(maybe_parse_next_tecplot_zone(reader, line_num));

      if (zone.has_value())
      {
        zones.emplace_back(zone.value());
        n_ignored_line += line_num - temp_line_num - 1 - zone.value().nnode();
      }
      else
      {
        std::cout << "done" << std::endl;
        break;
      }
    }

    auto nbzones = zones.size();
    std::cout << "* " + std::to_string(n_ignored_line) + " lines ignored"
              << std::endl;
    if (nbzones == 0)
    {
      //---------------------------------------------------------------------------
      // Dummy output
      //---------------------------------------------------------------------------
      std::ofstream writer(output_filename);
      if (!writer.is_open())
      {
        throw std::runtime_error("Failed to open output file for writing!");
      }

      writer.close();
      std::runtime_error(
        "The input file " + input_filename.string() +
        "\ndoes not contain any Tecplot zone! Created an empty file..."
        "\nYou may want to try converting this file to point data"
        "\nwith -p2 option if dim == 2 or -p3 option if dim == 3");
    }

    //---------------------------------------------------------------------------
    // Compute global informations
    //---------------------------------------------------------------------------
    // Get the solution dimension (compute the maximum value)
    uint8_t dimension = 0;
    for (const auto& zone : zones)
    {
      if (zone.Dimension[0] > dimension)
      {
        dimension = zone.Dimension[0];
      }
    }

    if (dimension == 1)
    {
      std::runtime_error(
        "1D is not supported. Use GnuPlot to display the data file.");
    }

    // Loop over the zones to get the number of nodes and cells
    uint64_t n_node = 0;
    uint64_t n_cell = 0;
    for (const auto& zone : zones)
    {
      n_node += zone.nnode();
      n_cell += zone.N_cell[0];
    }

    // Get the field count
    // Assumption: the number of fields is assumed to be constant over all
    // zones
    auto n_field = zones[0].Nodes[0].Fields.size();

    //---------------------------------------------------------------------------
    // Write into Vtk XML output file
    //---------------------------------------------------------------------------
    std::ofstream writer(output_filename);
    if (!writer.is_open())
    {
      throw std::runtime_error("Failed to open output file for writing!");
    }

    writer << VtkXml::Header << std::endl;
    writer << VtkXml::UnstructuredGridHeader << std::endl;
    writer << VtkXml::PieceHeader(n_node, n_cell) << std::endl;

    //---------------------------------------------------------------------------
    // Nodes
    //---------------------------------------------------------------------------
    std::cout << "Write nodal coordinates..................." << std::flush;
    writer << VtpXml::PointsHeader << std::endl;
    for (const auto& zone : zones)
    {
      for (const auto& node : zone.Nodes)
      {
        writer << std::scientific << node.coordinates[0] << " "
               << node.coordinates[1] << " " << node.coordinates[2] << " "
               << std::defaultfloat << std::endl;
      }
    }
    writer << VtpXml::PointsFooter << std::endl;
    std::cout << "done" << std::endl;

    //---------------------------------------------------------------------------
    // Cells
    //---------------------------------------------------------------------------
    writer << VtkXml::CellsHeader << std::endl;

    {
      //---------------------------------------------------------------------------
      // Cell connectivity
      //---------------------------------------------------------------------------
      std::cout << "Write nodal connectivity.................." << std::flush;
      writer << VtkXml::ConnectivityHeader << std::endl;

      uint64_t pos = 0; // Current cell origin node index
      uint64_t n_zone = 0;
      for (auto& zone : zones)
      {
        n_zone += 1;
        uint64_t max_node = 0;
        if (zone.Cell_format[0] == ZoneFormat::TriangleOrTetrahedron)
        {
          if (n_zone == 1)
          {
            max_node = 0;
          }

          // Dump connectivities
          for (uint64_t e = 0; e < zone.N_cell[0]; e++)
          {
            Connectivity& conn = zone.Connectivities[e];
            auto conn_size = conn.size();
            for (unsigned i = 0; i < conn_size; i++)
            {
              // renumber connectivity table if it belongs to zone
              // > 1
              // TODO: Overload in-place operators
              if (n_zone > 1)
              {
                conn[i] = conn[i] + max_node + 1;
              }
              conn[i] = conn[i] - 1;
              writer << conn[i] << " " << std::endl;
            }
            writer << std::endl;
          }

          // Compute maximum node number for renumbering of next
          // connectivity table
          max_node = 0;
          for (uint64_t e = 0; e < zone.N_cell[0]; e++)
          {
            auto& conn = zone.Connectivities[e];
            auto conn_count = conn.size();
            auto dummy =
              *std::max_element(conn.data().begin(), conn.data().end());
            if (dummy > max_node)
            {
              max_node = dummy;
            }
          }
        }

        auto write_vector =
          [](std::ofstream& writer, const std::vector<uint64_t>& indices)
        {
          auto n_index = indices.size();
          for (uint64_t i = 0; i < n_index - 1; i++)
          {
            writer << indices[i] << " ";
          }
          writer << indices[n_index - 1] << std::endl;
        };

        if (zone.Cell_format[0] == ZoneFormat::Other)
        {
          auto dimI = zone.Edges[0];
          std::vector<uint64_t> indices = {pos, pos + 1};

          // Line (dim I = 2)
          if (zone.Dimension[0] == 1)
          {
            indices = {pos, pos + 1};
            write_vector(writer, indices);
            pos += dimI;
          }

          // Quad
          if (zone.Dimension[0] == 2)
          {
            auto dimJ = zone.Edges[1];
            std::vector<uint64_t> indices(4, 0);
            for (uint64_t j = 0; j < dimJ - 1; j++)
            {
              for (uint64_t i = 0; i < dimI - 1; i++)
              {
                // Unique face of the cell
                indices[0] = pos; // bottom-left node
                indices[1] = pos + 1; // bottom-right node
                indices[2] = pos + 1 + dimI; // top-right node
                indices[3] = pos + dimI; // to-left node
                write_vector(writer, indices);
                // Next cell
                pos += 1;
              }
              // Next row of cells
              pos += 1;
            }
            // Next zone
            pos += dimI;
          }

          // Hexahedron
          if (zone.Dimension[0] == 3)
          {
            auto dimJ = zone.Edges[1];
            auto dimK = zone.Edges[2];
            std::vector<uint64_t> indices(8, 0);
            for (uint64_t k = 0; k < dimK - 1; k++)
            {
              for (uint64_t j = 0; j < dimJ - 1; j++)
              {
                for (uint64_t i = 0; i < dimI - 1; i++)
                {
                  // Front face of the cell
                  indices[0] = pos; // bottom-left node
                  indices[1] = pos + 1; // bottom-right node
                  indices[2] = pos + 1 + dimI; // top-right node
                  indices[3] = pos + dimI; // to-left node
                  // Back face of the cell
                  auto back_pos = pos + dimI * dimJ;
                  indices[4] = back_pos;
                  indices[5] = back_pos + 1;
                  indices[6] = back_pos + 1 + dimI;
                  indices[7] = back_pos + dimI;
                  write_vector(writer, indices);
                  // Next cell
                  pos += 1;
                }
                // Next row of cells
                pos += 1;
              }
              // Next k
              pos += dimI;
            }
            // Next zone
            pos += dimI * dimJ;
          }
        }
      }
      writer << VtkXml::ConnectivityFooter << std::endl;
      std::cout << "done" << std::endl;

      //---------------------------------------------------------------------------
      // Cell offset
      //---------------------------------------------------------------------------
      std::cout << "Write nodal offsets......................." << std::flush;
      writer << VtkXml::OffsetsHeader << std::endl;
      uint64_t offset = 0;
      for (const auto& zone : zones)
      {
        for (unsigned i = 0; i < zone.N_cell[0]; i++)
        {
          if (zone.Cell_type[0] == VTKCellType::Line) offset += 2;
          if (zone.Cell_type[0] == VTKCellType::Triangle) offset += 3;
          if (zone.Cell_type[0] == VTKCellType::Quad) offset += 4;
          if (zone.Cell_type[0] == VTKCellType::Tetrahedron) offset += 4;
          if (zone.Cell_type[0] == VTKCellType::Hexahedron) offset += 8;
          writer << offset << std::endl;
        }
      }
      writer << VtkXml::OffsetsFooter << std::endl;
      std::cout << "done" << std::endl;

      //---------------------------------------------------------------------------
      // Cell types
      //---------------------------------------------------------------------------
      std::cout << "Write nodal types........................." << std::flush;
      writer << VtkXml::TypesHeader << std::endl;
      auto cell_type0 = zones[0].Cell_type[0];
      bool warn = false;
      for (const auto& zone : zones)
      {
        std::string cell_type_as_str = std::to_string(zone.Cell_type[0]);
        if ((!warn) && (zone.Cell_type[0] != cell_type0))
        {
          warn = true;
        }
        writer << StringHelpers::repeat(cell_type_as_str + "\n", zone.N_cell[0])
               << std::flush;
      }
      writer << VtkXml::TypesFooter << std::endl;
      std::cout << "done" << std::endl;
      if (warn)
      {
        std::cout << "Warning: Different types of elements" << std::endl;
      }
    }

    writer << VtkXml::CellsFooter << std::endl;

    //---------------------------------------------------------------------------
    // Fields
    //---------------------------------------------------------------------------
    writer << VtkXml::PointDataHeader << std::endl;
    for (unsigned i_field = 0; i_field < n_field; i_field++)
    {
      printf("Write field %02d/%02lu.........................",
             i_field + 1,
             n_field);
      std::cout << std::flush;
      writer << VtkXml::FieldHeader(i_field + 1) << std::endl;
      for (const auto& zone : zones)
      {
        for (const auto& node : zone.Nodes)
        {
          writer << std::scientific << node.Fields[i_field] << std::defaultfloat
                 << std::endl;
        }
      }
      std::cout << "done" << std::endl;
      writer << VtkXml::FieldFooter << std::endl;
    }
    writer << VtkXml::PointDataFooter << std::endl;
    writer << VtkXml::PieceFooter << std::endl;
    writer << VtkXml::UnstructuredGridFooter << std::endl;
    writer << VtkXml::Footer << std::endl;
  }

  /************************************************************************************
   * @brief Converts Tecplot file generates by oomph-lib into Vtp XML file.
   *
   * @param input_filename:
   * @param output_filename:
   * @param dim:
   ************************************************************************************/
  void tecplot_to_vtpxml(const fs::path& input_filename,
                         const fs::path& output_filename,
                         const unsigned& dim)
  {
    //---------------------------------------------------------------------------
    // Retrieve the points from the input Tecplot file
    //---------------------------------------------------------------------------
    std::ifstream reader(input_filename);
    if (!reader.is_open())
    {
      throw std::runtime_error("Failed to open input file for reading!");
    }

    std::cout << "Parse input file for points........" << std::endl;
    uint64_t line_num{0};
    uint64_t prev_line_num{0};
    std::vector<int64_t> offset_list{};
    int64_t count{0};
    uint64_t n_zone{1};
    std::vector<TecplotNode> points{};

    while (true)
    {
      auto point =
        std::move(maybe_parse_next_tecplot_node(reader, line_num, dim));
      if (line_num != prev_line_num + 1)
      {
        if (prev_line_num != 0)
        {
          offset_list.emplace_back(count);
          n_zone += 1;
        }
      }
      count += 1;
      prev_line_num = line_num;

      if (point.has_value())
      {
        points.emplace_back(point.value());
      }
      else
      {
        std::cout << "done" << std::endl;
        break;
      }
    }

    offset_list.emplace_back(count - 1);
    auto nbpoints = points.size();

    if (nbpoints == 0)
    {
      throw std::runtime_error("The input file does not contain any points!");
    }

    //---------------------------------------------------------------------------
    // Compute global informations
    //---------------------------------------------------------------------------

    // Get the field count
    // Assumption: the number of fields is assumed to be constant over all
    // points
    auto n_field = points[0].Fields.size();

    //---------------------------------------------------------------------------
    // Write into Vtp XML output file
    //---------------------------------------------------------------------------
    std::ofstream writer(output_filename);
    if (!writer.is_open())
    {
      throw std::runtime_error("Failed to open output file for writing!");
    }

    writer << VtpXml::Header << std::endl;
    writer << VtpXml::PolyDataHeader << std::endl;
    writer << VtpXml::PieceHeader(nbpoints, n_zone) << std::endl;

    //---------------------------------------------------------------------------
    // Nodes
    //---------------------------------------------------------------------------
    std::cout << "Write points coordinates...................";
    writer << VtpXml::PointsHeader << std::endl;
    for (const auto& point : points)
    {
      writer << std::scientific << point.coordinates[0] << " "
             << point.coordinates[1] << " " << point.coordinates[2] << " "
             << std::defaultfloat << std::endl;
    }
    writer << VtpXml::PointsFooter << std::endl;
    std::cout << "done" << std::endl;

    //---------------------------------------------------------------------------
    // Fields
    //---------------------------------------------------------------------------
    writer << VtpXml::PointDataHeader << std::endl;
    for (unsigned i_field = 0; i_field < n_field; i_field++)
    {
      printf("Write field %02d/%02lu.........................",
             i_field + 1,
             n_field);
      std::cout << std::flush;
      writer << VtpXml::FieldHeader(i_field + 1) << std::endl;
      for (const auto& point : points)
      {
        writer << std::scientific << point.Fields[i_field] << std::defaultfloat
               << std::endl;
      }
      std::cout << "done" << std::endl;
      writer << VtpXml::FieldFooter << std::endl;
    }
    writer << VtpXml::PointDataFooter << std::endl;

    //---------------------------------------------------------------------------
    // Headers and footers
    //---------------------------------------------------------------------------

    writer << VtpXml::VertsHeader << std::endl;
    writer << VtpXml::VertsFooter << std::endl;
    writer << VtpXml::LinesHeader << std::endl;
    // Prepare line information:
    writer << VtkXml::ConnectivityHeader << std::endl;
    for (uint64_t i = 0; i < nbpoints; i++)
    {
      writer << i << " ";
    }
    writer << VtkXml::ConnectivityFooter << std::endl;
    writer << VtkXml::OffsetsHeader << std::endl;
    for (const auto offset : offset_list)
    {
      writer << offset << " ";
    }
    writer << VtkXml::OffsetsFooter << std::endl;
    // end line information
    writer << VtpXml::LinesFooter << std::endl;
    writer << VtpXml::StripsHeader << std::endl;
    writer << VtpXml::StripsFooter << std::endl;
    writer << VtpXml::PolysHeader << std::endl;
    writer << VtpXml::PolysFooter << std::endl;
    writer << VtpXml::PieceFooter << std::endl;
    writer << VtpXml::PolyDataFooter << std::endl;
    writer << VtpXml::Footer << std::endl;
  }
} // namespace Tecplot


/************************************************************************************
 * @brief Check if the output file exists, if so decide if we should overwrite
 *        it and tell the user what we are doing.
 *
 * @return true:
 * @return false:
 ************************************************************************************/
bool is_ok_to_write_to_output_file(const bool& flag,
                                   const fs::path& in_file_name,
                                   const fs::path& out_file_name)
{
  if (!fs::exists(out_file_name))
  {
    return true;
  }

  std::cout << "File " << out_file_name << " already exists!" << std::endl;

  // Get modification times (in seconds since unix epoch I think)
  auto in_mtime = std::filesystem::last_write_time(in_file_name);
  auto out_mtime = std::filesystem::last_write_time(out_file_name);

  if (flag)
  {
    std::cout << "Overwriting regardless of modification times because of flag."
              << std::endl;
    return true;
  }
  else if (in_mtime > out_mtime)
  {
    std::cout << "Overwriting because input file is newer than output file"
              << std::endl;
    return true;
  }
  else
  {
    std::cout << "Not overwriting." << std::endl;
  }
  return false;
}

/************************************************************************************
 * @brief
 *
 * @return int:
 ************************************************************************************/
int main(int argc, char** argv)
{
  std::cout << "* oomph-convert (C++), ver. 20221120" << std::endl;
  std::string input_file;
  std::string output_file;
  unsigned dim = -1;
  bool write_points_flag = false;
  bool overwrite_output_file = false;
  bool add_zero_padding_to_output_filename = false;

  using namespace oomph;

  // Parse commandline flags
  CommandLineArgs::setup(argc, argv);
  CommandLineArgs::specify_command_line_flag("-i", &input_file);
  CommandLineArgs::specify_command_line_flag("-o", &output_file);
  CommandLineArgs::specify_command_line_flag("-p", &dim);
  CommandLineArgs::specify_command_line_flag("-f");
  CommandLineArgs::specify_command_line_flag("-h");
  CommandLineArgs::specify_command_line_flag("-z");
  CommandLineArgs::parse_and_assign();

  // Print help and finish here
  if (CommandLineArgs::command_line_flag_has_been_set("-h"))
  {
    print_usage();
    exit(0);
  }

  // Make sure the input file exists
  if (input_file.empty())
  {
    print_usage();
    exit(0);
  }
  else if (!fs::exists(input_file))
  {
    throw std::runtime_error("File '" + input_file + "' does not exist!");
  }

  // Check the user asked for 2D or 3D output points
  if (dim != -1)
  {
    if ((dim != 2) && (dim != 3))
    {
      print_usage();
      exit(0);
    }
    bool write_points_flag = true;
  }
  if (CommandLineArgs::command_line_flag_has_been_set("-z"))
  {
    add_zero_padding_to_output_filename = true;
  }
  if (CommandLineArgs::command_line_flag_has_been_set("-f"))
  {
    overwrite_output_file = true;
  }

  // If the user didn't specify an output file
  if (output_file.empty())
  {
    // Construct the output filename ourselves
    if (write_points_flag)
    {
      output_file = fs::path(input_file).replace_extension(".vtp");
    }
    else
    {
      output_file = fs::path(input_file).replace_extension(".vtu");
    }
    std::cout << "Made own file: " << output_file << std::endl;
  }
  // If the user did specify an output file, lets make sure we can write to it
  else
  {
    if (fs::exists(output_file) && (!overwrite_output_file))
    {
      throw std::runtime_error("ERROR: Output File '" + output_file +
                               "' already exists!");
    }
  }

  // Get suffix of input/output files
  std::string isuffix = fs::path(input_file).extension();
  std::string osuffix = fs::path(output_file).extension();

  // Zero pad output names if requested
  // FIXME: not implemented yet
  if (add_zero_padding_to_output_filename)
  {
    auto add_trailing_zeros =
      [](const std::string& output_file,
         const std::string& output_suffix) -> std::string
    {
      throw std::runtime_error("Not implemented yet!");
      return "";
    };
    output_file = add_trailing_zeros(output_file, osuffix);
  }

  if (is_ok_to_write_to_output_file(
        overwrite_output_file, input_file, output_file))
  {
    using std::chrono::duration;
    using std::chrono::duration_cast;
    using std::chrono::high_resolution_clock;
    using std::chrono::milliseconds;

    auto t_start = high_resolution_clock::now();

    // Convert from oomph-lib Tecplot format to VTK XML format
    if ((isuffix == ".dat") and (osuffix == ".vtu"))
    {
      Tecplot::tecplot_to_vtkxml(input_file, output_file);
    }
    // Convert from oomph-lib Tecplot format to VTP XML format
    else if ((isuffix == ".dat") and (osuffix == ".vtp"))
    {
      Tecplot::tecplot_to_vtpxml(input_file, output_file, dim);
    }
    else
    {
      std::runtime_error("Sorry, cannot convert between " + isuffix + " and " +
                         osuffix + " file formats.");
    }
    auto t_end = high_resolution_clock::now();

    // Get number of seconds as a double
    auto t_ms = duration_cast<milliseconds>(t_end - t_start);
    double t_secs = t_ms.count() / 1000.0;

    std::cout << "* Conversion done in " << t_secs << " seconds" << std::endl;
    std::cout << "* Output file name: " << output_file << std::endl;
  }
}