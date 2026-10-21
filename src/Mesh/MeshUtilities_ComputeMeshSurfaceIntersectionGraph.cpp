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
#include <stack>

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
          std::map<unsigned int, unsigned int> signed_vertices;
          std::map<unsigned int, unsigned int> vertices;
          std::vector<unsigned int> edges;
          std::vector<std::list<unsigned int>> vertices_adjacency;
      };

      bool is_only_negative;
      bool is_only_positive;
      bool is_only_on_surface;

      graph negative_graph;
      graph positive_graph;
  };
  inline signed_graphs create_signed_graphs(const std::vector<unsigned int>& selected_vertices,
                                            const std::vector<unsigned int>& selected_edges,
                                            const std::vector<Gedim::MeshUtilities::MeshSurfaceIntersectionGraphData::VertexData>& graph_vertices,
                                            const std::vector<Gedim::MeshUtilities::MeshSurfaceIntersectionGraphData::EdgeData>& graph_edges)
  {
    signed_graphs result;

    result.is_only_negative = true;
    result.is_only_positive = true;
    result.is_only_on_surface = true;

    for (const auto& v : selected_vertices)
    {
      switch (graph_vertices.at(v).surface_type)
      {
        case Gedim::MeshUtilities::PointSurfacePositions::Positive:
        {
          result.positive_graph.vertices.insert(std::make_pair(v, result.positive_graph.vertices.size()));
          result.positive_graph.signed_vertices.insert(std::make_pair(v, result.positive_graph.signed_vertices.size()));
          result.is_only_negative = false;
          result.is_only_on_surface = false;
        }
          break;
        case Gedim::MeshUtilities::PointSurfacePositions::Negative:
        {
          result.negative_graph.vertices.insert(std::make_pair(v, result.negative_graph.vertices.size()));
          result.negative_graph.signed_vertices.insert(std::make_pair(v, result.negative_graph.signed_vertices.size()));
          result.is_only_positive = false;
          result.is_only_on_surface = false;
        }
          break;
        case Gedim::MeshUtilities::PointSurfacePositions::OnSurface:
        {
          result.positive_graph.vertices.insert(std::make_pair(v, result.positive_graph.vertices.size()));
          result.negative_graph.vertices.insert(std::make_pair(v, result.negative_graph.vertices.size()));
        }
          continue;
        default:
          throw std::runtime_error("Unknown position type");
      }
    }

    result.negative_graph.edges.reserve(selected_edges.size());
    result.positive_graph.edges.reserve(selected_edges.size());

    result.negative_graph.vertices_adjacency.resize(result.negative_graph.vertices.size());
    result.positive_graph.vertices_adjacency.resize(result.positive_graph.vertices.size());

    for (const auto& e : selected_edges)
    {
      const auto& graph_edge = graph_edges.at(e);
      const unsigned int origin = graph_edge.extremes.at(0);
      const unsigned int end = graph_edge.extremes.at(1);

      switch (graph_edge.surface_type)
      {
        case Gedim::MeshUtilities::PointSurfacePositions::Positive:
        {
          result.positive_graph.edges.push_back(e);
          const unsigned int origin_v = result.positive_graph.vertices.at(origin);
          const unsigned int end_v = result.positive_graph.vertices.at(end);

          result.positive_graph.vertices_adjacency.at(origin_v).push_back(end);
          result.positive_graph.vertices_adjacency.at(end_v).push_back(origin);
        }
          break;
        case Gedim::MeshUtilities::PointSurfacePositions::Negative:
        {
          result.negative_graph.edges.push_back(e);
          const unsigned int origin_v = result.negative_graph.vertices.at(origin);
          const unsigned int end_v = result.negative_graph.vertices.at(end);

          result.negative_graph.vertices_adjacency.at(origin_v).push_back(end);
          result.negative_graph.vertices_adjacency.at(end_v).push_back(origin);
        }
          break;
        case Gedim::MeshUtilities::PointSurfacePositions::OnSurface:
        {
          result.positive_graph.edges.push_back(e);
          const unsigned int pos_origin_v = result.positive_graph.vertices.at(origin);
          const unsigned int pos_end_v = result.positive_graph.vertices.at(end);

          result.positive_graph.vertices_adjacency.at(pos_origin_v).push_back(end);
          result.positive_graph.vertices_adjacency.at(pos_end_v).push_back(origin);

          result.negative_graph.edges.push_back(e);
          const unsigned int neg_origin_v = result.negative_graph.vertices.at(origin);
          const unsigned int neg_end_v = result.negative_graph.vertices.at(end);

          result.negative_graph.vertices_adjacency.at(neg_origin_v).push_back(end);
          result.negative_graph.vertices_adjacency.at(neg_end_v).push_back(origin);
        }
          break;
        default:
          throw std::runtime_error("Unknown position type");
      }
    }

    result.negative_graph.edges.shrink_to_fit();
    result.positive_graph.edges.shrink_to_fit();

    return result;
  }
  // ***************************************************************************
  struct new_face_data final
  {
      std::vector<unsigned int> vertices;
      std::vector<unsigned int> edges;
  };
  inline std::vector<new_face_data> create_faces_from_graph(const signed_graphs::graph& graph,
                                                            const std::vector<Gedim::MeshUtilities::MeshSurfaceIntersectionGraphData::VertexData>& graph_vertices,
                                                            const std::vector<Gedim::MeshUtilities::MeshSurfaceIntersectionGraphData::EdgeData>& graph_edges)
  {
    if (graph_vertices.size() < 3 ||
        graph_edges.size() < 3)
      return {};

    std::vector<new_face_data> result;

    result.reserve(graph_vertices.size() - 3);

    std::vector<bool> vertices_included_in_face(graph.signed_vertices.size(),
                                                false);
    for (const auto& vertex : graph.signed_vertices)
    {
      std::vector<bool> vertices_visited(graph.vertices.size(),
                                         false);

      for (const auto& s_v : graph.signed_vertices)
      {
        const auto s_v_idx = graph.vertices.at(s_v.first);
        vertices_visited.at(s_v_idx) = vertices_included_in_face.at(s_v.second);
      }

      std::vector<int> vertices_prev(graph.vertices.size(),
                                     -1);
      std::stack<unsigned int> stack;

      int prev_vert = -1;
      stack.push(vertex.first);
      while (!stack.empty())
      {
        const auto visited = stack.top();
        const auto visited_idx = graph.vertices.at(visited);
        vertices_visited.at(visited_idx) = true;
        vertices_prev.at(visited_idx) = prev_vert;
        stack.pop();

        for (const auto adjacent : graph.vertices_adjacency[visited_idx])
        {
          if (adjacent == vertex.first &&
              vertices_prev.at(visited_idx) != static_cast<int>(vertex.first))
          {
            vertices_prev.at(vertex.second) = visited;
            break;
          }

          if (vertices_visited.at(graph.vertices.at(adjacent)))
            continue;

          stack.push(adjacent);
        }

        prev_vert = visited;

        if (vertices_prev.at(vertex.second) != -1)
          break;
      }

      if (vertices_prev.at(vertex.second) == -1)
        throw std::runtime_error("New faces not found for vertex " + std::to_string(vertex.first));

      std::list<unsigned int> new_face_vertices;
      unsigned int visited = vertex.first;
      do
      {
        new_face_vertices.push_back(visited);
        visited = vertices_prev.at(graph.vertices.at(visited));
      }
      while (vertices_prev.at(graph.vertices.at(visited)) !=
             static_cast<int>(vertex.first));
      new_face_vertices.push_back(visited);

      result.push_back({});
      auto& new_face = result.back();
      new_face.vertices.reserve(new_face_vertices.size());
      for(auto iter = new_face_vertices.rbegin();
          iter != new_face_vertices.rend();
          iter++)
      {
        new_face.vertices.push_back(*iter);

        const auto signed_vertex = graph.signed_vertices.find(*iter);
        if (signed_vertex != graph.signed_vertices.end())
          vertices_included_in_face.at(signed_vertex->second) = true;
      }
    }

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
      const auto cell2D_edges = mesh.Cell2DEdges(f);

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

      for (const auto& s : surface_segments.cell2Ds_surface_segments.at(f))
        face_edges.push_back(result.surface_segment_to_edge.at(s));

      face_vertices.shrink_to_fit();
      face_edges.shrink_to_fit();

      const auto face_graph = create_signed_graphs(face_vertices,
                                                   face_edges,
                                                   result.vertices,
                                                   result.edges);

      if (face_graph.is_only_negative ||
          face_graph.is_only_positive ||
          face_graph.is_only_on_surface)
      {
        const auto face_sign = face_graph.is_only_on_surface ?
                                 Gedim::MeshUtilities::PointSurfacePositions::OnSurface :
                                 face_graph.is_only_negative ?
                                   Gedim::MeshUtilities::PointSurfacePositions::Negative :
                                   Gedim::MeshUtilities::PointSurfacePositions::Positive;

        result.mesh_cell2D_to_faces.at(f) = result.faces.size();
        result.faces.push_back(
              {
                face_vertices,
                face_edges,
                face_sign,
                f,
                Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::None
              });
        continue;
      }

      const auto new_negative_faces = create_faces_from_graph(face_graph.negative_graph,
                                                              result.vertices,
                                                              result.edges);
      const auto new_positive_faces = create_faces_from_graph(face_graph.positive_graph,
                                                              result.vertices,
                                                              result.edges);
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
