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
  inline unsigned int find_edge(const std::array<unsigned int, 2>& edge_extremes,
                                const std::map<std::array<unsigned int, 2>, unsigned int>& edges_map)
  {
    std::array<unsigned int, 2> map_extremes =
    {
      std::min(edge_extremes.at(0), edge_extremes.at(1)),
      std::max(edge_extremes.at(0), edge_extremes.at(1))
    };

    const auto edge_found = edges_map.find(map_extremes);

    if (edge_found == edges_map.end())
      throw std::runtime_error("graph edge not found!");

    return edge_found->second;
  }
  // ***************************************************************************
  inline std::vector<unsigned int> find_face_edges(const std::vector<unsigned int>& face_vertices,
                                                   const std::map<std::array<unsigned int, 2>, unsigned int>& edges_map)
  {
    std::vector<unsigned int> face_edges(face_vertices.size());
    for (unsigned int v = 0; v < face_vertices.size(); ++v)
    {
      face_edges.at(v) = find_edge({
                                     face_vertices.at(v),
                                     face_vertices.at((v + 1) % face_vertices.size())
                                   },
                                   edges_map);
    }

    return face_edges;
  }
  // ***************************************************************************
  inline unsigned int create_edge(const std::array<unsigned int, 2>& edge_extremes,
                                  const Gedim::MeshUtilities::PointSurfacePositions& edge_surface_position,
                                  const unsigned int mesh_cell1D_index,
                                  const unsigned int surface_segment_index,
                                  Gedim::MeshUtilities::MeshSurfaceIntersectionGraphData& result)
  {
    unsigned int edge_index = result.edges.size();

    std::array<unsigned int, 2> map_extremes =
    {
      std::min(edge_extremes.at(0), edge_extremes.at(1)),
      std::max(edge_extremes.at(0), edge_extremes.at(1))
    };

    result.edges_map.insert(std::make_pair(map_extremes,
                                           edge_index));
    result.edges.push_back(
          {
            edge_extremes,
            edge_surface_position,
            mesh_cell1D_index,
            surface_segment_index
          });

    return edge_index;
  }
  // ***************************************************************************
  struct signed_graphs final
  {
      struct graph final
      {
          std::vector<unsigned int> vertices;
          std::vector<unsigned int> edges;
      };

      graph negative_graph;
      graph positive_graph;
  };
  signed_graphs create_signed_graphs(const std::vector<unsigned int>& selected_vertices,
                                     const std::vector<unsigned int>& selected_edges,
                                     const std::vector<Gedim::MeshUtilities::MeshSurfaceIntersectionGraphData::VertexData>& graph_vertices,
                                     const std::vector<Gedim::MeshUtilities::MeshSurfaceIntersectionGraphData::EdgeData>& graph_edges)
  {
    signed_graphs result;
    result.negative_graph.vertices.reserve(selected_vertices.size());
    result.positive_graph.vertices.reserve(selected_vertices.size());
    result.negative_graph.edges.reserve(selected_edges.size());
    result.positive_graph.edges.reserve(selected_edges.size());

    for (const auto& v : selected_vertices)
    {
      switch (graph_vertices.at(v).surface_type)
      {
        case Gedim::MeshUtilities::PointSurfacePositions::Positive:
          result.positive_graph.vertices.push_back(v);
          break;
        case Gedim::MeshUtilities::PointSurfacePositions::Negative:
          result.negative_graph.vertices.push_back(v);
          break;
        case Gedim::MeshUtilities::PointSurfacePositions::OnSurface:
          result.positive_graph.vertices.push_back(v);
          result.negative_graph.vertices.push_back(v);
          continue;
        default:
          throw std::runtime_error("Unknown position type");
      }
    }

    for (const auto& e : selected_edges)
    {
      switch (graph_edges.at(e).surface_type)
      {
        case Gedim::MeshUtilities::PointSurfacePositions::Positive:
          result.positive_graph.edges.push_back(e);
          break;
        case Gedim::MeshUtilities::PointSurfacePositions::Negative:
          result.negative_graph.edges.push_back(e);
          break;
        case Gedim::MeshUtilities::PointSurfacePositions::OnSurface:
          result.positive_graph.edges.push_back(e);
          result.negative_graph.edges.push_back(e);
          continue;
        default:
          throw std::runtime_error("Unknown position type");
      }
    }

    result.negative_graph.vertices.shrink_to_fit();
    result.positive_graph.vertices.shrink_to_fit();
    result.negative_graph.edges.shrink_to_fit();
    result.positive_graph.edges.shrink_to_fit();

    return result;
  }
  // ***************************************************************************
  Gedim::MeshUtilities::MeshSurfaceIntersectionGraphData Gedim::MeshUtilities::ComputeMeshSurfaceIntersectionGraph(
      const Gedim::MeshUtilities::IntersectMeshWithSurfaceResult &surface_intersections,
      const Gedim::MeshUtilities::MarchingSquaresResult::SurfaceSegmentsData &surface_segments,
      const Gedim::MeshUtilities::MarchingCubesResult::SurfaceTrianglesData& surface_triangles,
      const Gedim::IMeshDAO &mesh)
  {
    Gedim::MeshUtilities::MeshSurfaceIntersectionGraphData result;

    result.mesh_cell0D_to_vertex.resize(mesh.Cell0DTotalNumber());
    result.surface_point_to_vertex.resize(surface_intersections.surface_intersection_points.size());
    result.vertices.reserve(mesh.Cell0DTotalNumber() + surface_intersections.surface_intersection_points.size());

    for (unsigned int v = 0; v < mesh.Cell0DTotalNumber(); ++v)
    {
      const auto vertex_surface_point = surface_intersections.cell0Ds_surface_point.at(v);

      switch (vertex_surface_point.surface_position)
      {
        case Gedim::MeshUtilities::PointSurfacePositions::Positive:
        case Gedim::MeshUtilities::PointSurfacePositions::Negative: {
          result.mesh_cell0D_to_vertex.at(v) = result.vertices.size();
          result.vertices.push_back({mesh.Cell0DCoordinates(v),
                                     surface_intersections.cell0Ds_surface_point.at(v).surface_position,
                                     Gedim::MeshUtilities::MeshSurfaceIntersectionGraphData::VertexData::MeshTypes::Vertex,
                                     v,
                                     surface_intersections.cell0Ds_surface_point.at(v).index});
        }
          break;
        case Gedim::MeshUtilities::PointSurfacePositions::OnSurface:
          continue;
        default:
          throw std::runtime_error("Unknown position type");
      }
    }

    for (unsigned int s = 0; s < surface_intersections.surface_intersection_points.size(); ++s)
    {
      const auto &surface_intersection = surface_intersections.surface_intersection_points.at(s);

      switch (surface_intersection.Type)
      {
        case Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::SurfacePointData::Types::Vertex: {
          result.mesh_cell0D_to_vertex.at(surface_intersection.mesh_id) = result.vertices.size();
          result.surface_point_to_vertex.at(s) = result.vertices.size();
          result.vertices.push_back({surface_intersection.coordinate,
                                     Gedim::MeshUtilities::PointSurfacePositions::OnSurface,
                                     Gedim::MeshUtilities::MeshSurfaceIntersectionGraphData::VertexData::MeshTypes::Vertex,
                                     surface_intersection.mesh_id,
                                     s});
        }
          continue;
        case Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::SurfacePointData::Types::Edge: {
          result.surface_point_to_vertex.at(s) = result.vertices.size();
          result.vertices.push_back({surface_intersection.coordinate,
                                     Gedim::MeshUtilities::PointSurfacePositions::OnSurface,
                                     Gedim::MeshUtilities::MeshSurfaceIntersectionGraphData::VertexData::MeshTypes::Edge,
                                     surface_intersection.mesh_id,
                                     s});
        }
          break;
        default:
          throw std::runtime_error("Unknown intersection type");
      }
    }

    result.vertices.shrink_to_fit();

    result.mesh_cell1D_to_edges.resize(mesh.Cell1DTotalNumber());
    result.surface_segment_to_edge.resize(surface_segments.surface_segments.size());
    result.edges.reserve(2 * mesh.Cell1DTotalNumber() + surface_segments.surface_segments.size());

    for (unsigned int e = 0; e < mesh.Cell1DTotalNumber(); ++e)
    {
      const auto edge_origin_id = mesh.Cell1DOrigin(e);
      const auto edge_end_id = mesh.Cell1DEnd(e);

      const auto edge_surface_point = surface_intersections.cell1Ds_surface_point.at(e);

      switch (edge_surface_point.surface_position)
      {
        case Gedim::MeshUtilities::PointSurfacePositions::Positive:
        case Gedim::MeshUtilities::PointSurfacePositions::Negative:
        {
          const auto new_edge_id = create_edge({
                                                 result.mesh_cell0D_to_vertex.at(edge_origin_id),
                                                 result.mesh_cell0D_to_vertex.at(edge_end_id)
                                               },
                                               edge_surface_point.surface_position,
                                               e,
                                               Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::None,
                                               result);

          result.mesh_cell1D_to_edges.at(e) = { new_edge_id };
        }
          continue;
        case Gedim::MeshUtilities::PointSurfacePositions::OnSurface:
        {
          const auto edge_origin_surface_position = surface_intersections.cell0Ds_surface_point.at(edge_origin_id).surface_position;
          const auto edge_end_surface_position = surface_intersections.cell0Ds_surface_point.at(edge_end_id).surface_position;

          if (edge_origin_surface_position == Gedim::MeshUtilities::PointSurfacePositions::OnSurface &&
              edge_end_surface_position == Gedim::MeshUtilities::PointSurfacePositions::OnSurface)
            continue;

          const auto new_edge_id_one = create_edge({
                                                     result.mesh_cell0D_to_vertex.at(edge_origin_id),
                                                     result.surface_point_to_vertex.at(edge_surface_point.index),
                                                   },
                                                   edge_origin_surface_position,
                                                   e,
                                                   edge_surface_point.index,
                                                   result);

          const auto new_edge_id_two = create_edge({
                                                     result.surface_point_to_vertex.at(edge_surface_point.index),
                                                     result.mesh_cell0D_to_vertex.at(edge_end_id),
                                                   },
                                                   edge_end_surface_position,
                                                   e,
                                                   edge_surface_point.index,
                                                   result);

          result.mesh_cell1D_to_edges.at(e) =
          {
              new_edge_id_one,
              new_edge_id_two
          };

        }
          break;
        default:
          throw std::runtime_error("Unknown position type");
      }
    }

    for (unsigned int s = 0; s < surface_segments.surface_segments.size(); ++s)
    {
      const auto &segment_data = surface_segments.surface_segments.at(s);

      const auto new_edge_id = create_edge({
                                             result.surface_point_to_vertex.at(segment_data.indices.at(0)),
                                             result.surface_point_to_vertex.at(segment_data.indices.at(1))
                                           },
                                           Gedim::MeshUtilities::PointSurfacePositions::OnSurface,
                                           Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::None,
                                           s,
                                           result);

      result.surface_segment_to_edge.at(s) = new_edge_id;
    }

    result.edges.shrink_to_fit();

    result.mesh_cell2D_to_faces.resize(mesh.Cell2DTotalNumber());
    result.surface_triangle_to_face.resize(surface_triangles.surface_triangles.size());
    result.faces.reserve(3 * mesh.Cell2DTotalNumber() +
                         surface_triangles.surface_triangles.size());

    for (unsigned int f = 0; f < mesh.Cell2DTotalNumber(); ++f)
    {
      const auto cell2D_vertices = mesh.Cell2DVertices(f);
      const auto cell2D_edges = mesh.Cell2DVertices(f);

      std::vector<unsigned int> face_vertices;
      std::vector<unsigned int> face_edges;
      face_vertices.reserve(cell2D_vertices.size() +
                            cell2D_edges.size());
      face_edges.reserve(cell2D_vertices.size() +
                         cell2D_edges.size());

      for (const auto& cell0D_index : cell2D_vertices)
        face_vertices.push_back(result.mesh_cell0D_to_vertex.at(cell0D_index));

      for (const auto& cell1D_index : cell2D_edges)
      {
        const auto surface_point_index = surface_intersections.cell1Ds_surface_point.at(cell1D_index);

        if (surface_point_index.index !=
            Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::None)
          face_vertices.push_back(result.surface_point_to_vertex.at(surface_point_index.index));

        for (const auto& edge_index : result.mesh_cell1D_to_edges.at(cell1D_index))
          face_edges.push_back(edge_index);
      }

      face_vertices.shrink_to_fit();
      face_edges.shrink_to_fit();

      const auto face_graph = create_signed_graphs(face_vertices,
                                                   face_edges,
                                                   result.vertices,
                                                   result.edges);

      if (face_graph.negative_graph.vertices.empty())
      {
        result.mesh_cell2D_to_faces.at(f) = result.faces.size();
        result.faces.push_back(
              {
                face_vertices,
                face_edges,
                Gedim::MeshUtilities::PointSurfacePositions::Positive,
                f,
                Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::None
              });
        continue;
      }

      if (face_graph.positive_graph.vertices.empty())
      {
        result.mesh_cell2D_to_faces.at(f) = result.faces.size();
        result.faces.push_back(
              {
                face_vertices,
                face_edges,
                Gedim::MeshUtilities::PointSurfacePositions::Negative,
                f,
                Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::None
              });
        continue;
      }
    }

    for (unsigned int t = 0; t < surface_triangles.surface_triangles.size(); ++t)
    {
      const auto &triangle_data = surface_triangles.surface_triangles.at(t);

      const std::vector<unsigned int> face_vertices = {
        result.surface_point_to_vertex.at(triangle_data.indices.at(0)),
        result.surface_point_to_vertex.at(triangle_data.indices.at(1)),
        result.surface_point_to_vertex.at(triangle_data.indices.at(2))
      };
      const auto face_edges = find_face_edges(face_vertices,
                                              result.edges_map);

      result.surface_triangle_to_face.at(t) = result.faces.size();
      result.faces.push_back(
            {
              face_vertices,
              face_edges,
              Gedim::MeshUtilities::PointSurfacePositions::OnSurface,
              Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::None,
              t
            });
    }

    return result;
  }
  // ***************************************************************************
} // namespace Gedim
