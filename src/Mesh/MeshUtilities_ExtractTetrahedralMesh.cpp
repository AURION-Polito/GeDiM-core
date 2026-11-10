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
