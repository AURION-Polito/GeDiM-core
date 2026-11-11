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
  // **************************************************************edges*************
  constexpr std::array<std::array<unsigned int, 3>, 4> tetrahedron_faces =
  {
    std::array<unsigned int, 3>({ 0, 1, 2 }),
    std::array<unsigned int, 3>({ 0, 1, 3 }),
    std::array<unsigned int, 3>({ 0, 2, 3 }),
    std::array<unsigned int, 3>({ 1, 2, 3 })
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
  Gedim::MeshUtilities::ExtractTetrahedralMeshData Gedim::MeshUtilities::ExtractTriangularMesh(
      const Gedim::MeshUtilities::MeshCell2DsTrianglesData& cell2Ds_triangles,
      const Gedim::IMeshDAO &mesh)
  {
    return ExtractTetrahedralMesh(cell2Ds_triangles,
                                  {},
                                  mesh);
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
      result.cell2Ds_additional_point_to_vertex.at(f).shrink_to_fit();
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
      result.cell3Ds_additional_point_to_vertex.at(c).shrink_to_fit();
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

    std::map<std::array<unsigned int, 3>, unsigned int> faces;
    result.triangles.reserve(num_cell2D_triangles +
                             3 * num_cell3D_tetra);
    for (unsigned int f = 0; f < mesh.Cell2DTotalNumber(); ++f)
    {
      const auto& cell2D_triangles = cell2Ds_triangles.cell2Ds_triangles.at(f);
      result.original_cell2Ds_to_triangles.at(f).resize(cell2D_triangles.size());

      for (unsigned int t = 0; t < cell2D_triangles.size(); ++t)
      {
        const auto& cell2D_triangle = cell2D_triangles.at(t);

        Gedim::MeshUtilities::ExtractTetrahedralMeshData::TriangleData new_cell2D_triangle;
        new_cell2D_triangle.original_cell2D_index = f;

        for (unsigned int t_v = 0; t_v < 3; ++t_v)
        {
          const auto& local_v = cell2D_triangle.vertices.at(t_v);

          unsigned int vertex = Gedim::MeshUtilities::None;
          if (local_v.original_vertex_index != Gedim::MeshUtilities::None)
            vertex = result.original_cell0Ds_to_vertex.at(mesh.Cell2DVertex(f, local_v.original_vertex_index));
          else if (local_v.additional_vertex_index != Gedim::MeshUtilities::None)
            vertex = result.cell2Ds_additional_point_to_vertex.at(f).at(local_v.additional_vertex_index);
          else
            throw std::runtime_error("triangle vertex not recognized");

          new_cell2D_triangle.vertices.at(t_v) = vertex;
        }

        for (unsigned int t_v = 0; t_v < 3; ++t_v)
        {
          const auto& edge_o = new_cell2D_triangle.vertices.at(t_v);
          const auto& edge_e = new_cell2D_triangle.vertices.at((t_v + 1) % 3);

          std::array<unsigned int, 2> edge = {
            std::min(edge_o, edge_e),
            std::max(edge_o, edge_e)
          };

          const auto edge_found = edges.find(edge);

          if (edge_found == edges.end())
          {
            edges.insert(std::make_pair(edge, result.edges.size()));
            new_cell2D_triangle.edges.at(t_v) = result.edges.size();
            result.edges.push_back({
                                     edge,
                                     Gedim::MeshUtilities::None
                                   });
          }
          else
            new_cell2D_triangle.edges.at(t_v) = edge_found->second;
        }

        faces.insert(std::make_pair(new_cell2D_triangle.vertices, result.triangles.size()));
        result.original_cell2Ds_to_triangles.at(f).at(t) = result.triangles.size();
        result.triangles.push_back(new_cell2D_triangle);
      }
    }

    result.tetrahedra.reserve(num_cell3D_tetra);
    for (unsigned int c = 0; c < mesh.Cell3DTotalNumber(); ++c)
    {
      const auto& cell3D_tetrahedra = cell3Ds_tetrahedra.cell3Ds_tetrahedra.at(c);
      result.original_cell3Ds_to_tetrahedra.at(c).resize(cell3D_tetrahedra.size());

      for (unsigned int t = 0; t < cell3D_tetrahedra.size(); ++t)
      {
        Gedim::MeshUtilities::ExtractTetrahedralMeshData::TetrahedronData new_cell3D_tetrahedron;
        new_cell3D_tetrahedron.original_cell3D_index = c;
        const auto& cell3D_tetrahedron = cell3D_tetrahedra.at(t);

        for (unsigned int t_v = 0; t_v < 4; ++t_v)
        {
          const auto& local_v = cell3D_tetrahedron.vertices.at(t_v);

          unsigned int vertex = Gedim::MeshUtilities::None;
          if (local_v.original_vertex_index != Gedim::MeshUtilities::None)
            vertex = result.original_cell0Ds_to_vertex.at(mesh.Cell3DVertex(c, local_v.original_vertex_index));
          else if (local_v.additional_vertex_index != Gedim::MeshUtilities::None)
          {
            const auto& cell3D_additional_point = cell3Ds_tetrahedra.cell3Ds_additional_points.at(c).at(local_v.additional_vertex_index);

            if (cell3D_additional_point.face_index == Gedim::MeshUtilities::None)
              vertex = result.cell3Ds_additional_point_to_vertex.at(c).at(local_v.additional_vertex_index);
            else
            {
              const auto c_2D = mesh.Cell3DFace(c, cell3D_additional_point.face_index);
              vertex = result.cell2Ds_additional_point_to_vertex.at(c_2D).at(cell3D_additional_point.face_additional_point_index);
            }
          }
          else
            throw std::runtime_error("tetrahedron vertex not recognized");

          new_cell3D_tetrahedron.vertices.at(t_v) = vertex;
        }

        for (unsigned int e = 0; e < tetrahedron_edges.size(); ++e)
        {
          const auto& edge_o = new_cell3D_tetrahedron.vertices.at(tetrahedron_edges.at(e).at(0));
          const auto& edge_e = new_cell3D_tetrahedron.vertices.at(tetrahedron_edges.at(e).at(1));

          std::array<unsigned int, 2> edge = {
            std::min(edge_o, edge_e),
            std::max(edge_o, edge_e)
          };

          const auto edge_found = edges.find(edge);

          if (edge_found == edges.end())
          {
            edges.insert(std::make_pair(edge, result.edges.size()));
            new_cell3D_tetrahedron.edges.at(e) = result.edges.size();
            result.edges.push_back({
                                     edge,
                                     Gedim::MeshUtilities::None
                                   });
          }
          else
            new_cell3D_tetrahedron.edges.at(e) = edge_found->second;
        }

        for (unsigned int f = 0; f < tetrahedron_faces.size(); ++f)
        {
          Gedim::MeshUtilities::ExtractTetrahedralMeshData::TriangleData new_cell3D_triangle;

          new_cell3D_triangle.vertices = {
            new_cell3D_tetrahedron.vertices.at(tetrahedron_faces.at(f).at(0)),
            new_cell3D_tetrahedron.vertices.at(tetrahedron_faces.at(f).at(1)),
            new_cell3D_tetrahedron.vertices.at(tetrahedron_faces.at(f).at(2))
          };
          new_cell3D_triangle.original_cell2D_index = Gedim::MeshUtilities::None;

          const auto triangle_found = faces.find(new_cell3D_triangle.vertices);

          if (triangle_found == faces.end())
          {
            for (unsigned int t_v = 0; t_v < 3; ++t_v)
            {
              const auto& edge_o = new_cell3D_triangle.vertices.at(t_v);
              const auto& edge_e = new_cell3D_triangle.vertices.at((t_v + 1) % 3);

              std::array<unsigned int, 2> edge = {
                std::min(edge_o, edge_e),
                std::max(edge_o, edge_e)
              };

              new_cell3D_triangle.edges.at(t_v) = edges.at(edge);
            }

            faces.insert(std::make_pair(new_cell3D_triangle.vertices, result.triangles.size()));
            new_cell3D_tetrahedron.faces.at(f) = result.triangles.size();
            result.triangles.push_back(new_cell3D_triangle);
          }
          else
            new_cell3D_tetrahedron.faces.at(f) = triangle_found->second;
        }

        result.original_cell3Ds_to_tetrahedra.at(c).at(t) = result.tetrahedra.size();
        result.tetrahedra.push_back(new_cell3D_tetrahedron);
      }
    }

    result.edges.shrink_to_fit();
    result.triangles.shrink_to_fit();
    result.tetrahedra.shrink_to_fit();

    return result;
  }
  // ***************************************************************************
} // namespace Gedim
