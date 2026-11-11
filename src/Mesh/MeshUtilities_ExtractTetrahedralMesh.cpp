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
  constexpr std::array<std::array<unsigned int, 2>, 6> tetrahedron_edges =
  {
    std::array<unsigned int, 2>({ 0, 1 }),
    std::array<unsigned int, 2>({ 1, 2 }),
    std::array<unsigned int, 2>({ 0, 2 }),
    std::array<unsigned int, 2>({ 0, 3 }),
    std::array<unsigned int, 2>({ 1, 3 }),
    std::array<unsigned int, 2>({ 2, 3 })
  };
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
            triangles.at(t).vertices.at(0) = { triangulation_indices.at(3 * t), Gedim::MeshUtilities::None };
            triangles.at(t).vertices.at(1) = { triangulation_indices.at(3 * t + 1), Gedim::MeshUtilities::None };
            triangles.at(t).vertices.at(2) = { triangulation_indices.at(3 * t + 2), Gedim::MeshUtilities::None };
          };
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
          cell2D_additional_points = {
            {
              cell2D_barycenter
            }
          };

          for (unsigned int t = 0; t < numTriangles; t++)
          {
            triangles.at(t).vertices.at(0) = { Gedim::MeshUtilities::None, 0 };
            triangles.at(t).vertices.at(1) = { triangulation_indices.at(3 * t + 1), Gedim::MeshUtilities::None };
            triangles.at(t).vertices.at(2) = { triangulation_indices.at(3 * t + 2), Gedim::MeshUtilities::None };
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
          const unsigned int cell3D_barycenter_position = cell3D_additional_points.size();
          cell3D_additional_points.push_back({
                                               cell3D_barycenter,
                                               Gedim::MeshUtilities::None,
                                               Gedim::MeshUtilities::None,
                                             });
          for (unsigned int f = 0; f < mesh.Cell3DNumberFaces(c_3D); ++f)
          {
            const unsigned int c_2D = mesh.Cell3DFace(c_3D, f);
            const auto& cell2D_additional_points = cell2Ds_triangles.cell2Ds_additional_points.at(c_2D);
            const auto& cell2D_triangles = cell2Ds_triangles.cell2Ds_triangles.at(c_2D);

            for (unsigned int f_t = 0; f_t < cell2D_triangles.size(); ++f_t)
            {
              const auto& cell2D_triangle = cell2D_triangles.at(f_t);

              std::array<Gedim::MeshUtilities::MeshCell3DsTetrahedraData::TetrahedronData::Vertex, 4> tetrahedron;
              for (unsigned int f_t_v = 0; f_t_v < 3; ++f_t_v)
              {
                if (cell2D_triangle.vertices.at(f_t_v).additional_vertex_index != Gedim::MeshUtilities::None)
                {
                  const auto& cell2D_additional_point = cell2D_additional_points.at(cell2D_triangle.vertices.at(f_t_v).additional_vertex_index);

                  tetrahedron.at(f_t_v).original_vertex_index = Gedim::MeshUtilities::None;
                  tetrahedron.at(f_t_v).additional_vertex_index = cell3D_additional_points.size();
                  cell3D_additional_points.push_back({
                                                       cell2D_additional_point.coordinate,
                                                       f,
                                                       cell2D_triangle.vertices.at(f_t_v).additional_vertex_index
                                                     });
                }
                else if (cell2D_triangle.vertices.at(f_t_v).original_vertex_index != Gedim::MeshUtilities::None)
                {
                  const unsigned int cell0D_index = mesh.Cell2DVertex(c_2D, cell2D_triangle.vertices.at(f_t_v).original_vertex_index);
                  tetrahedron.at(f_t_v).original_vertex_index = mesh.Cell3DFindVertex(c_3D, cell0D_index);
                  tetrahedron.at(f_t_v).additional_vertex_index = Gedim::MeshUtilities::None;
                }
                else
                  throw std::runtime_error("Triangle point not recognized");
              }

              tetrahedron.at(3).original_vertex_index = Gedim::MeshUtilities::None;
              tetrahedron.at(3).additional_vertex_index = cell3D_barycenter_position;
              tetrahedra.push_back({ tetrahedron });

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
      const Gedim::MeshUtilities::MeshCell2DsTrianglesData& cell2Ds_triangles,
      const Gedim::MeshUtilities::MeshCell3DsTetrahedraData& cell3Ds_tetrahedra,
      const Gedim::IMeshDAO &mesh)
  {
    Gedim::MeshUtilities::ExtractTetrahedralMeshData result;

    result.original_cell0Ds_to_vertex.resize(mesh.Cell0DTotalNumber(),
                                             Gedim::MeshUtilities::None);
    result.original_cell1Ds_to_edge.resize(mesh.Cell1DTotalNumber(),
                                           Gedim::MeshUtilities::None);
    result.original_cell2Ds_to_triangles.resize(mesh.Cell2DTotalNumber());
    result.original_cell3Ds_to_tetrahedra.resize(mesh.Cell3DTotalNumber());

    result.cell2Ds_additional_point_to_vertex.resize(mesh.Cell2DTotalNumber());
    result.cell3Ds_additional_point_to_vertex.resize(mesh.Cell3DTotalNumber());

    result.vertices.reserve(mesh.Cell0DTotalNumber() +
                            cell2Ds_triangles.cell2Ds_additional_points.size() +
                            cell3Ds_tetrahedra.cell3Ds_additional_points.size());

    for (unsigned int v = 0; v < mesh.Cell0DTotalNumber(); ++v)
    {
      result.original_cell0Ds_to_vertex.at(v) = result.vertices.size();
      result.vertices.push_back({
                                  mesh.Cell0DCoordinates(v),
                                  v
                                });
    }

    unsigned int num_cell2D_triangles = 0;
    for (unsigned int f = 0; f < mesh.Cell2DTotalNumber(); ++f)
    {
      const auto& cell2D_triangles = cell2Ds_triangles.cell2Ds_triangles.at(f);
      num_cell2D_triangles += cell2D_triangles.size();

      result.cell2Ds_additional_point_to_vertex.at(f).reserve(cell2Ds_triangles.cell2Ds_additional_points.size());
      for (const auto& additional_point : cell2Ds_triangles.cell2Ds_additional_points.at(f))
      {
        result.cell2Ds_additional_point_to_vertex.at(f).push_back(result.vertices.size());
        result.vertices.push_back({
                                    additional_point.coordinate,
                                    Gedim::MeshUtilities::None
                                  });
      }
    }

    unsigned int num_cell3D_tetra = 0;
    for (unsigned int c = 0; c < mesh.Cell3DTotalNumber(); ++c)
    {
      const auto& cell3D_tetra = cell3Ds_tetrahedra.cell3Ds_tetrahedra.at(c);
      num_cell3D_tetra += cell3D_tetra.size();

      result.cell3Ds_additional_point_to_vertex.at(c).reserve(cell3Ds_tetrahedra.cell3Ds_additional_points.size());
      for (const auto& additional_point : cell3Ds_tetrahedra.cell3Ds_additional_points.at(c))
      {
        if (additional_point.face_index != Gedim::MeshUtilities::None)
          continue;

        result.cell3Ds_additional_point_to_vertex.at(c).push_back(result.vertices.size());
        result.vertices.push_back({
                                    additional_point.coordinate,
                                    Gedim::MeshUtilities::None
                                  });
      }
    }
    result.vertices.shrink_to_fit();

    std::map<std::array<unsigned int, 2>, unsigned int> edges;
    result.edges.reserve(3 * num_cell2D_triangles +
                         3 * num_cell3D_tetra);

    for (unsigned int e = 0; e < mesh.Cell1DTotalNumber(); ++e)
    {
      result.original_cell1Ds_to_edge.at(e) = result.edges.size();

      const auto& extremes = mesh.Cell1DExtremes(e);

      std::array<unsigned int, 2> edge = {
        static_cast<unsigned int>(std::min(extremes[0], extremes[1])),
        static_cast<unsigned int>(std::max(extremes[0], extremes[1]))
      };

      edges.insert(std::make_pair(edge, result.edges.size()));
      result.edges.push_back({
                               edge,
                               e
                             });
    }

    for (unsigned int f = 0; f < mesh.Cell2DTotalNumber(); ++f)
    {
      const auto& cell2D_triangles = cell2Ds_triangles.cell2Ds_triangles.at(f);
      for (unsigned int t = 0; t < cell2D_triangles.size(); ++t)
      {
        const auto& cell2D_triangle = cell2D_triangles.at(t);

        for (unsigned int t_v = 0; t_v < 3; ++t_v)
        {
          const auto& local_edge_o = cell2D_triangle.vertices.at(t_v);
          const auto& local_edge_e = cell2D_triangle.vertices.at((t_v + 1) % 3);

          unsigned int edge_o = Gedim::MeshUtilities::None;
          if (local_edge_o.original_vertex_index != Gedim::MeshUtilities::None)
            edge_o = result.original_cell0Ds_to_vertex.at(mesh.Cell2DVertex(f, local_edge_o.original_vertex_index));
          else if (local_edge_o.additional_vertex_index != Gedim::MeshUtilities::None)
            edge_o = result.cell2Ds_additional_point_to_vertex.at(f).at(local_edge_o.additional_vertex_index);
          else
            throw std::runtime_error("triangle vertex not recognized");

          unsigned int edge_e = Gedim::MeshUtilities::None;
          if (local_edge_e.original_vertex_index != Gedim::MeshUtilities::None)
            edge_e = result.original_cell0Ds_to_vertex.at(mesh.Cell2DVertex(f, local_edge_e.original_vertex_index));
          else if (local_edge_e.additional_vertex_index != Gedim::MeshUtilities::None)
            edge_e = result.cell2Ds_additional_point_to_vertex.at(f).at(local_edge_e.additional_vertex_index);
          else
            throw std::runtime_error("triangle vertex not recognized");

          std::array<unsigned int, 2> edge = {
            std::min(edge_o, edge_e),
            std::max(edge_o, edge_e)
          };

          if (!edges.contains(edge))
          {
            edges.insert(std::make_pair(edge, result.edges.size()));
            result.edges.push_back({
                                     edge,
                                     Gedim::MeshUtilities::None
                                   });
          }
        }
      }
    }

    for (unsigned int c = 0; c < mesh.Cell3DTotalNumber(); ++c)
    {
      const auto& cell3D_tetrahedra = cell3Ds_tetrahedra.cell3Ds_tetrahedra.at(c);
      for (unsigned int t = 0; t < cell3D_tetrahedra.size(); ++t)
      {
        const auto& cell3D_tetrahedron = cell3D_tetrahedra.at(t);

        for (unsigned int e = 0; e < tetrahedron_edges.size(); ++e)
        {
          const auto& local_edge_o = cell3D_tetrahedron.vertices.at(tetrahedron_edges.at(e).at(0));
          const auto& local_edge_e = cell3D_tetrahedron.vertices.at(tetrahedron_edges.at(e).at(1));

          unsigned int edge_o = Gedim::MeshUtilities::None;
          if (local_edge_o.original_vertex_index != Gedim::MeshUtilities::None)
            edge_o = result.original_cell0Ds_to_vertex.at(mesh.Cell3DVertex(c, local_edge_o.original_vertex_index));
          else if (local_edge_o.additional_vertex_index != Gedim::MeshUtilities::None)
          {
            const auto& cell3D_additional_point = cell3Ds_tetrahedra.cell3Ds_additional_points.at(c).at(local_edge_o.additional_vertex_index);

            if (cell3D_additional_point.face_index == Gedim::MeshUtilities::None)
              edge_o = result.cell3Ds_additional_point_to_vertex.at(c).at(local_edge_o.additional_vertex_index);
            else
            {
              const auto c_2D = mesh.Cell3DFace(c, cell3D_additional_point.face_index);
              edge_o = result.cell2Ds_additional_point_to_vertex.at(c_2D).at(cell3D_additional_point.face_additional_point_index);
            }
          }
          else
            throw std::runtime_error("tetrahedron vertex not recognized");


          unsigned int edge_e = Gedim::MeshUtilities::None;
          if (local_edge_e.original_vertex_index != Gedim::MeshUtilities::None)
            edge_e = result.original_cell0Ds_to_vertex.at(mesh.Cell3DVertex(c, local_edge_e.original_vertex_index));
          else if (local_edge_e.additional_vertex_index != Gedim::MeshUtilities::None)
          {
            const auto& cell3D_additional_point = cell3Ds_tetrahedra.cell3Ds_additional_points.at(c).at(local_edge_e.additional_vertex_index);

            if (cell3D_additional_point.face_index == Gedim::MeshUtilities::None)
              edge_e = result.cell3Ds_additional_point_to_vertex.at(c).at(local_edge_e.additional_vertex_index);
            else
            {
              const auto c_2D = mesh.Cell3DFace(c, cell3D_additional_point.face_index);
              edge_e = result.cell2Ds_additional_point_to_vertex.at(c_2D).at(cell3D_additional_point.face_additional_point_index);
            }
          }
          else
            throw std::runtime_error("tetrahedron vertex not recognized");

          std::array<unsigned int, 2> edge = {
            std::min(edge_o, edge_e),
            std::max(edge_o, edge_e)
          };

          if (!edges.contains(edge))
          {
            edges.insert(std::make_pair(edge, result.edges.size()));
            result.edges.push_back({
                                     edge,
                                     Gedim::MeshUtilities::None
                                   });
          }
        }
      }
    }
    result.edges.shrink_to_fit();

    result.original_cell2Ds_to_triangles.reserve(num_cell2D_triangles +
                                                 3 * num_cell3D_tetra);
    result.original_cell3Ds_to_tetrahedra.reserve(num_cell3D_tetra);

    return result;
  }
  // ***************************************************************************
} // namespace Gedim
