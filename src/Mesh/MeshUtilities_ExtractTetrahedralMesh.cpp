// _LICENSE_HEADER_
//
// Copyright (C) 2019 - 2025.
// Terms register on the GPL-3.0 license.
//
// This file can be redistributed and/or modified under the license terms.
//
// See top level LICENSE file for more details.
//
// This file can be used citing references in CITATION.cff file.

#include "MeshUtilities.hpp"

namespace Gedim
{
  // ***************************************************************************
  Gedim::MeshUtilities::MeshCell2DsTrianglesData Gedim::MeshUtilities::ComputeMeshCell2DsTriangles(const Gedim::GeometryUtilities& geometry_utilities,
                                                                                                   const Gedim::MeshUtilities::MeshCell2DsTrianglesData::TrianglesTypes triangles_type,
                                                                                                   const Gedim::IMeshDAO &mesh)
  {
    Gedim::MeshUtilities::MeshCell2DsTrianglesData result;

    result.cell2Ds_additional_points.resize(mesh.Cell2DTotalNumber());
    result.cell2Ds_triangles.resize(mesh.Cell2DTotalNumber());

    switch (triangles_type)
    {
      case Gedim::MeshUtilities::MeshCell2DsTrianglesData::TrianglesTypes::ByFirstVertex:
      {
        for (unsigned int c_2D = 0; c_2D < mesh.Cell2DTotalNumber(); ++c_2D)
        {
          auto& triangles = result.cell2Ds_triangles.at(c_2D);

          const auto cell2D_vertices = mesh.Cell2DVerticesCoordinates(c_2D);
          const auto triangulation_indices = geometry_utilities.PolygonTriangulationByFirstVertex(cell2D_vertices);

          const unsigned int numTriangles = triangulation_indices.size() / 3;
          triangles.resize(numTriangles);

          for (unsigned int t = 0; t < numTriangles; t++)
          {
            triangles.at(t) = {
                              triangulation_indices.at(3 * t),
                              triangulation_indices.at(3 * t + 1),
                              triangulation_indices.at(3 * t + 2) };
          }
        }
      }
        break;
      case Gedim::MeshUtilities::MeshCell2DsTrianglesData::TrianglesTypes::ByBarycenter:
      {
        for (unsigned int c_2D = 0; c_2D < mesh.Cell2DTotalNumber(); ++c_2D)
        {
          auto& cell2D_additional_points = result.cell2Ds_additional_points.at(c_2D);
          auto& triangles = result.cell2Ds_triangles.at(c_2D);

          const auto cell2D_vertices = mesh.Cell2DVerticesCoordinates(c_2D);
          const auto cell2D_barycenter = geometry_utilities.PolygonBarycenter(cell2D_vertices);
          const auto triangulation_indices = geometry_utilities.PolygonTriangulationByInternalPoint(cell2D_vertices,
                                                                                                    cell2D_barycenter);

          const unsigned int numTriangles = triangulation_indices.size() / 3;
          triangles.resize(numTriangles);
          cell2D_additional_points = { cell2D_barycenter };

          for (unsigned int t = 0; t < numTriangles; t++)
          {
            triangles.at(t) = {
                              mesh.Cell2DNumberVertices(c_2D),
                              triangulation_indices.at(3 * t + 1),
                              triangulation_indices.at(3 * t + 2) };
          }
        }
      }
        break;
      default:
        throw std::runtime_error("Unknown triangles type");
    }

    return result;
  }
  // ***************************************************************************
  Gedim::MeshUtilities::MeshCell3DsTetrahedraData Gedim::MeshUtilities::ComputeMeshCell3DsTetrahedra(const Gedim::GeometryUtilities& geometry_utilities,
                                                                                                     const Gedim::MeshUtilities::MeshCell2DsTrianglesData& cell2Ds_triangles,
                                                                                                     const Gedim::MeshUtilities::MeshCell3DsTetrahedraData::TetrahedraTypes tetrahedra_type,
                                                                                                     const Gedim::IMeshDAO &mesh)
  {
    Gedim::MeshUtilities::MeshCell3DsTetrahedraData result;


    result.cell3Ds_additional_points.resize(mesh.Cell3DTotalNumber());
    result.cell3Ds_tetrahedra.resize(mesh.Cell3DTotalNumber());

    switch (tetrahedra_type)
    {
      case Gedim::MeshUtilities::MeshCell3DsTetrahedraData::TetrahedraTypes::ByBarycenter:
      {
        for (unsigned int c_3D = 0; c_3D < mesh.Cell3DTotalNumber(); ++c_3D)
        {
          auto& cell3D_additional_points = result.cell3Ds_additional_points.at(c_3D);
          auto& tetrahedra = result.cell3Ds_tetrahedra.at(c_3D);

          const auto cell3D_vertices = mesh.Cell3DVerticesCoordinates(c_3D);
          const auto cell3D_barycenter = geometry_utilities.PolyhedronBarycenter(cell3D_vertices);

          cell3D_additional_points.reserve(mesh.Cell3DNumberFaces(c_3D) + 1);

          unsigned int num_tetrahedra = 0;
          for (unsigned int f = 0; f < mesh.Cell3DNumberFaces(c_3D); ++f)
            num_tetrahedra += cell2Ds_triangles.cell2Ds_triangles.at(mesh.Cell3DFace(c_3D, f)).size();

          tetrahedra.reserve(num_tetrahedra);
          cell3D_additional_points.push_back(cell3D_barycenter);
          const unsigned int cell3D_barycenter_position = mesh.Cell3DNumberVertices(c_3D) + cell3D_additional_points.size();
          for (unsigned int f = 0; f < mesh.Cell3DNumberFaces(c_3D); ++f)
          {
            const unsigned int c_2D = mesh.Cell3DFace(c_3D, f);
            const auto& cell2D_additional_points = cell2Ds_triangles.cell2Ds_additional_points.at(c_2D);
            const auto& cell2D_triangles = cell2Ds_triangles.cell2Ds_triangles.at(c_2D);

            for (const auto& cell2D_triangle : cell2D_triangles)
            {
              std::array<unsigned int, 4> tetrahedron;
              for (unsigned int f_t_v = 0; f_t_v < 3; ++f_t_v)
              {
                if (cell2D_triangle.at(f_t_v) >= mesh.Cell2DNumberVertices(c_2D))
                {
                  tetrahedron.at(f_t_v) = mesh.Cell3DNumberVertices(c_3D) + cell3D_additional_points.size();
                  cell3D_additional_points.push_back(cell2D_additional_points.at(cell2D_triangle.at(f_t_v) - mesh.Cell2DNumberVertices(c_2D)));
                }
                else
                {
                  const unsigned int cell0D_index = mesh.Cell2DVertex(c_2D, cell2D_triangle.at(f_t_v));
                  tetrahedron.at(f_t_v) = mesh.Cell3DFindVertex(c_3D, cell0D_index);
                }
              }

              tetrahedron.at(3) = cell3D_barycenter_position;
              tetrahedra.push_back(tetrahedron);

            }
          }
          cell3D_additional_points.shrink_to_fit();
          tetrahedra.shrink_to_fit();
        }
      }
        break;
      default:
        throw std::runtime_error("Unknown tetrahedra type");
    }

    return result;
  }
  // ***************************************************************************
  Gedim::MeshUtilities::ExtractTetrahedralMeshData Gedim::MeshUtilities::ExtractTetrahedralMesh(
      const std::vector<std::vector<std::array<unsigned int, 3>>> &cell2Ds_triangles,
      const std::vector<std::vector<Eigen::Matrix3d>> &cell2Ds_triangles_3D_vertices,
      const std::vector<std::vector<std::array<unsigned int, 4>>> &cell3Ds_tetrahedra,
      const std::vector<std::vector<Eigen::MatrixXd>> &cell3Ds_tetrahedra_vertices,
      const Gedim::IMeshDAO &mesh)
  {
    Gedim::MeshUtilities::ExtractTetrahedralMeshData result;

    result.original_cell0Ds_to_vertex.resize(mesh.Cell0DTotalNumber(),
                                             Gedim::MeshUtilities::None);
    result.original_cell1Ds_to_edge.resize(mesh.Cell1DTotalNumber(),
                                           Gedim::MeshUtilities::None);
    result.original_cell2Ds_to_triangles.resize(mesh.Cell2DTotalNumber());
    result.original_cell3Ds_to_tetrahedra.resize(mesh.Cell3DTotalNumber());

    result.vertices.reserve(mesh.Cell0DTotalNumber() +
                            cell2Ds_triangles.size() +
                            cell3Ds_tetrahedra.size());

    for (unsigned int v = 0; v < mesh.Cell0DTotalNumber(); ++v)
    {
      result.original_cell0Ds_to_vertex.at(v) = result.vertices.size();
      result.vertices.push_back({
                                  mesh.Cell0DCoordinates(v),
                                  v
                                });
    }

    unsigned int max_num_cell2D_triangles = 0;
    for (unsigned int f = 0; f < mesh.Cell2DTotalNumber(); ++f)
    {
      const unsigned int cell2D_num_vertices = mesh.Cell2DNumberVertices(f);

      const auto& cell2D_triangles = cell2Ds_triangles.at(f);
      if (max_num_cell2D_triangles < cell2D_triangles.size())
        max_num_cell2D_triangles = cell2D_triangles.size();

      for (unsigned int t = 0; t < cell2D_triangles.size(); ++t)
      {
        const auto& triangle = cell2D_triangles.at(t);
        for (unsigned int t_v = 0; t_v < 3; ++t_v)
        {
          if (triangle.at(t_v) > cell2D_num_vertices)
          {
            result.vertices.push_back({
                                        cell2Ds_triangles_3D_vertices.at(f).at(t).col(t_v),
                                        Gedim::MeshUtilities::None
                                      });
          }
        }
      }
    }

    unsigned int max_num_cell3D_tetra = 0;
    for (unsigned int c = 0; c < mesh.Cell3DTotalNumber(); ++c)
    {
      const unsigned int cell3D_num_vertices = mesh.Cell3DNumberVertices(c);

      const auto& cell3D_tetra = cell3Ds_tetrahedra.at(c);
      if (max_num_cell3D_tetra < cell3D_tetra.size())
        max_num_cell3D_tetra = cell3D_tetra.size();

      for (unsigned int t = 0; t < cell3D_tetra.size(); ++t)
      {
        const auto& tetra = cell3D_tetra.at(t);
        for (unsigned int t_v = 0; t_v < 4; ++t_v)
        {
          if (tetra.at(t_v) > cell3D_num_vertices)
          {
            result.vertices.push_back({
                                        cell3Ds_tetrahedra_vertices.at(c).at(t).col(t_v),
                                        Gedim::MeshUtilities::None
                                      });
          }
        }
      }
    }

    result.edges.reserve(mesh.Cell1DTotalNumber() +
                         max_num_cell2D_triangles * cell2Ds_triangles.size() +
                         max_num_cell3D_tetra * cell3Ds_tetrahedra.size());
    result.original_cell2Ds_to_triangles.reserve(max_num_cell2D_triangles *
                                                 mesh.Cell2DTotalNumber());
    result.original_cell3Ds_to_tetrahedra.reserve(max_num_cell3D_tetra *
                                                  mesh.Cell3DTotalNumber());

    return result;
  }
  // ***************************************************************************
} // namespace Gedim
