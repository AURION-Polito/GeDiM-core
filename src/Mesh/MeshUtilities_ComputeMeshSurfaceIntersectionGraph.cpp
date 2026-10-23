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
#include <queue>
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
          std::map<unsigned int, unsigned int> signed_faces;
          std::map<unsigned int, unsigned int> faces;

          struct adj_vertex_element final
          {
              unsigned int vertex;
              unsigned int edge;
          };

          std::vector<std::list<adj_vertex_element>> vertices_edges_adjacency;

          struct adj_face_element final
          {
              unsigned int face;
              unsigned int edge;
          };

          std::vector<std::list<adj_face_element>> faces_edges_adjacency;
      };

      bool is_only_negative;
      bool is_only_positive;
      bool is_only_on_surface;

      graph negative_graph;
      graph positive_graph;
  };
  inline signed_graphs create_signed_graphs(const std::vector<unsigned int>& selected_vertices,
                                            const std::map<unsigned int, unsigned int>& selected_edges,
                                            const std::map<unsigned int, unsigned int>& selected_faces,
                                            const std::vector<Gedim::MeshUtilities::MeshSurfaceIntersectionGraphData::VertexData>& graph_vertices,
                                            const std::vector<Gedim::MeshUtilities::MeshSurfaceIntersectionGraphData::EdgeData>& graph_edges,
                                            const std::vector<Gedim::MeshUtilities::MeshSurfaceIntersectionGraphData::FaceData>& graph_faces)
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

    result.negative_graph.vertices_edges_adjacency.resize(result.negative_graph.vertices.size());
    result.positive_graph.vertices_edges_adjacency.resize(result.positive_graph.vertices.size());

    for (const auto& e : selected_edges)
    {
      const auto& graph_edge = graph_edges.at(e.first);
      const unsigned int origin = graph_edge.extremes.at(0);
      const unsigned int end = graph_edge.extremes.at(1);

      switch (graph_edge.surface_type)
      {
        case Gedim::MeshUtilities::PointSurfacePositions::Positive:
        {
          result.positive_graph.edges.push_back(e.first);
          const unsigned int origin_v = result.positive_graph.vertices.at(origin);
          const unsigned int end_v = result.positive_graph.vertices.at(end);

          result.positive_graph.vertices_edges_adjacency.at(origin_v).push_back({
                                                                                  end,
                                                                                  e.first
                                                                                });
          result.positive_graph.vertices_edges_adjacency.at(end_v).push_back({
                                                                               origin,
                                                                               e.first
                                                                             });
        }
          break;
        case Gedim::MeshUtilities::PointSurfacePositions::Negative:
        {
          result.negative_graph.edges.push_back(e.first);
          const unsigned int origin_v = result.negative_graph.vertices.at(origin);
          const unsigned int end_v = result.negative_graph.vertices.at(end);

          result.negative_graph.vertices_edges_adjacency.at(origin_v).push_back({
                                                                                  end,
                                                                                  e.first
                                                                                });
          result.negative_graph.vertices_edges_adjacency.at(end_v).push_back({
                                                                               origin,
                                                                               e.first
                                                                             });
        }
          break;
        case Gedim::MeshUtilities::PointSurfacePositions::OnSurface:
        {
          result.positive_graph.edges.push_back(e.first);
          const unsigned int pos_origin_v = result.positive_graph.vertices.at(origin);
          const unsigned int pos_end_v = result.positive_graph.vertices.at(end);

          result.positive_graph.vertices_edges_adjacency.at(pos_origin_v).push_back({
                                                                                      end,
                                                                                      e.first
                                                                                    });
          result.positive_graph.vertices_edges_adjacency.at(pos_end_v).push_back({
                                                                                   origin,
                                                                                   e.first
                                                                                 });

          result.negative_graph.edges.push_back(e.first);
          const unsigned int neg_origin_v = result.negative_graph.vertices.at(origin);
          const unsigned int neg_end_v = result.negative_graph.vertices.at(end);

          result.negative_graph.vertices_edges_adjacency.at(neg_origin_v).push_back({
                                                                                      end,
                                                                                      e.first
                                                                                    });
          result.negative_graph.vertices_edges_adjacency.at(neg_end_v).push_back({
                                                                                   origin,
                                                                                   e.first
                                                                                 });
        }
          break;
        default:
          throw std::runtime_error("Unknown position type");
      }
    }

    result.negative_graph.edges.shrink_to_fit();
    result.positive_graph.edges.shrink_to_fit();

    if (selected_faces.size() > 1)
    {

      std::vector<std::vector<unsigned int>> edges_faces(selected_edges.size());
      for (unsigned int e = 0; e < selected_edges.size(); ++e)
        edges_faces.at(e).reserve(selected_faces.size());


      for (const auto& f : selected_faces)
      {
        const auto& graph_face = graph_faces.at(f.first);

        for (const auto f_e : graph_face.intersection_graph_edges)
        {
          edges_faces.at(selected_edges.at(f_e)).push_back(f.first);
        }

        switch (graph_face.surface_type)
        {
          case Gedim::MeshUtilities::PointSurfacePositions::Positive:
          {
            result.positive_graph.faces.insert(std::make_pair(f.first, result.positive_graph.faces.size()));
            result.positive_graph.signed_faces.insert(std::make_pair(f.first, result.positive_graph.signed_faces.size()));
          }
            break;
          case Gedim::MeshUtilities::PointSurfacePositions::Negative:
          {
            result.negative_graph.faces.insert(std::make_pair(f.first, result.negative_graph.faces.size()));
            result.negative_graph.signed_faces.insert(std::make_pair(f.first, result.negative_graph.signed_faces.size()));
          }
            break;
          case Gedim::MeshUtilities::PointSurfacePositions::OnSurface:
          {
            result.positive_graph.faces.insert(std::make_pair(f.first, result.positive_graph.faces.size()));
            result.negative_graph.faces.insert(std::make_pair(f.first, result.negative_graph.faces.size()));
          }
            break;
          default:
            throw std::runtime_error("Unknown position type");
        }
      }

      for (unsigned int e = 0; e < selected_edges.size(); ++e)
        edges_faces.at(e).shrink_to_fit();

      result.negative_graph.faces_edges_adjacency.resize(result.negative_graph.faces.size());
      result.positive_graph.faces_edges_adjacency.resize(result.positive_graph.faces.size());

      for (const auto& e : selected_edges)
      {
        for (unsigned int e_f = 0; e_f < edges_faces.at(e.second).size(); ++e_f)
        {
          const auto f_one = edges_faces.at(e.second).at(e_f);
          const auto f_two = edges_faces.at(e.second).at((e_f + 1) % edges_faces.at(e.second).size());

          const auto& graph_face_one = graph_faces.at(f_one);
          const auto& graph_face_two = graph_faces.at(f_two);

          if (graph_face_one.surface_type == Gedim::MeshUtilities::PointSurfacePositions::OnSurface ||
              graph_face_one.surface_type == Gedim::MeshUtilities::PointSurfacePositions::Positive)
          {
            if (graph_face_two.surface_type == Gedim::MeshUtilities::PointSurfacePositions::OnSurface ||
                graph_face_two.surface_type == Gedim::MeshUtilities::PointSurfacePositions::Positive)
            {
              const auto f_one_idx = result.positive_graph.faces.at(f_one);
              const auto f_two_idx = result.positive_graph.faces.at(f_two);

              result.positive_graph.faces_edges_adjacency.at(f_one_idx).push_back({
                                                                                    f_two,
                                                                                    e.first
                                                                                  });
              result.positive_graph.faces_edges_adjacency.at(f_two_idx).push_back({
                                                                                    f_one,
                                                                                    e.first
                                                                                  });
            }
          }

          if (graph_face_one.surface_type == Gedim::MeshUtilities::PointSurfacePositions::OnSurface ||
              graph_face_one.surface_type == Gedim::MeshUtilities::PointSurfacePositions::Negative)
          {
            if (graph_face_two.surface_type == Gedim::MeshUtilities::PointSurfacePositions::OnSurface ||
                graph_face_two.surface_type == Gedim::MeshUtilities::PointSurfacePositions::Negative)
            {
              const auto f_one_idx = result.negative_graph.faces.at(f_one);
              const auto f_two_idx = result.negative_graph.faces.at(f_two);

              result.negative_graph.faces_edges_adjacency.at(f_one_idx).push_back({
                                                                                    f_two,
                                                                                    e.first
                                                                                  });
              result.negative_graph.faces_edges_adjacency.at(f_two_idx).push_back({
                                                                                    f_one,
                                                                                    e.first
                                                                                  });
            }
          }
        }
      }
    }

    return result;
  }
  // ***************************************************************************
  struct new_face_data final
  {
      std::vector<unsigned int> vertices;
      std::vector<unsigned int> edges;
  };
  inline Gedim::GeometryUtilities::PolygonOrientations check_new_face_orientation(const std::map<unsigned int, unsigned int>& original_face_vertex_to_mesh_cell1D,
                                                                                  const unsigned int new_graph_vertex,
                                                                                  const unsigned int new_graph_edge,
                                                                                  const std::vector<Gedim::MeshUtilities::MeshSurfaceIntersectionGraphData::EdgeData>& graph_edges)
  {
    const auto& original_vertex = original_face_vertex_to_mesh_cell1D.find(new_graph_vertex);
    if (original_vertex != original_face_vertex_to_mesh_cell1D.end())
    {
      const auto& new_edge = graph_edges.at(new_graph_edge);
      if (new_edge.original_mesh_cell1D_index != Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::None)
      {
        if (new_edge.original_mesh_cell1D_index == original_vertex->second)
          return Gedim::GeometryUtilities::PolygonOrientations::Clockwise;
        else
          return Gedim::GeometryUtilities::PolygonOrientations::CounterClockwise;
      }
    }

    return Gedim::GeometryUtilities::PolygonOrientations::Unknown;
  }
  inline std::vector<new_face_data> create_faces_from_graph(const std::map<unsigned int, unsigned int>& original_face_vertex_to_mesh_cell1D,
                                                            const signed_graphs::graph& graph,
                                                            const std::vector<Gedim::MeshUtilities::MeshSurfaceIntersectionGraphData::EdgeData>& graph_edges)
  {
    if (graph.vertices.size() < 3 ||
        graph.edges.size() < 3)
      return {};

    std::vector<new_face_data> result;

    result.reserve(graph.vertices.size());

    std::vector<bool> vertices_included_in_face(graph.signed_vertices.size(),
                                                false);

    for (const auto& vertex : graph.signed_vertices)
    {
      if (vertices_included_in_face.at(vertex.second))
        continue;

      const unsigned int vertex_idx = graph.vertices.at(vertex.first);

      std::vector<bool> vertices_visited(graph.vertices.size(),
                                         false);

      for (const auto& s_v : graph.signed_vertices)
      {
        const auto s_v_idx = graph.vertices.at(s_v.first);
        vertices_visited.at(s_v_idx) = vertices_included_in_face.at(s_v.second);
      }

      std::vector<std::array<int, 2>> vertices_prev(graph.vertices.size(),
                                                    { -1, -1 });
      std::stack<std::array<int, 2>> stack;

      std::array<int, 2> prev_vert = { -1, -1 };
      stack.push({ static_cast<int>(vertex.first), -1 });
      while (!stack.empty())
      {
        const auto visited = stack.top();
        const auto visited_idx = graph.vertices.at(visited.at(0));
        vertices_visited.at(visited_idx) = true;
        vertices_prev.at(visited_idx) = { prev_vert.at(0), visited.at(1) };
        stack.pop();

        for (const auto adjacent : graph.vertices_edges_adjacency[visited_idx])
        {
          if (adjacent.vertex == vertex.first &&
              vertices_prev.at(visited_idx).at(0) != static_cast<int>(vertex.first))
          {
            vertices_prev.at(vertex_idx) = { visited.at(0), static_cast<int>(adjacent.edge) };
            break;
          }

          if (vertices_visited.at(graph.vertices.at(adjacent.vertex)))
            continue;

          stack.push({
                       static_cast<int>(adjacent.vertex),
                       static_cast<int>(adjacent.edge)
                     });
        }

        prev_vert = visited;

        if (vertices_prev.at(vertex_idx).at(0) != -1)
          break;
      }

      if (vertices_prev.at(vertex_idx).at(0) == -1 ||
          vertices_prev.at(vertex_idx).at(1) == -1)
        throw std::runtime_error("New faces not found for vertex " + std::to_string(vertex.first));


      Gedim::GeometryUtilities::PolygonOrientations new_face_direction = Gedim::GeometryUtilities::PolygonOrientations::Unknown;
      std::list<std::array<int, 2>> new_face_list;
      std::array<int, 2> visited = vertices_prev.at(graph.vertices.at(vertex.first));
      do
      {
        new_face_list.push_back(visited);
        if (new_face_direction == Gedim::GeometryUtilities::PolygonOrientations::Unknown)
          new_face_direction = check_new_face_orientation(original_face_vertex_to_mesh_cell1D,
                                                          visited.at(0),
                                                          visited.at(1),
                                                          graph_edges);


        visited = vertices_prev.at(graph.vertices.at(visited.at(0)));
      }
      while (visited.at(0) !=
             static_cast<int>(vertex.first));
      new_face_list.push_back(visited);
      if (new_face_direction == Gedim::GeometryUtilities::PolygonOrientations::Unknown)
        new_face_direction = check_new_face_orientation(original_face_vertex_to_mesh_cell1D,
                                                        visited.at(0),
                                                        visited.at(1),
                                                        graph_edges);

      if (new_face_direction == Gedim::GeometryUtilities::PolygonOrientations::Unknown)
        throw std::runtime_error("New faces direction not found for vertex " + std::to_string(vertex.first));

      result.push_back({});
      auto& new_face = result.back();
      new_face.vertices.resize(new_face_list.size());
      new_face.edges.resize(new_face_list.size());

      if (new_face_direction == Gedim::GeometryUtilities::PolygonOrientations::Clockwise)
      {
        unsigned int fv = 0;
        for(auto iter = new_face_list.rbegin();
            iter != new_face_list.rend();
            iter++)
        {
          new_face.vertices.at(fv) = iter->at(0);
          new_face.edges.at(fv) = iter->at(1);

          const auto signed_vertex = graph.signed_vertices.find(iter->at(0));
          if (signed_vertex != graph.signed_vertices.end())
            vertices_included_in_face.at(signed_vertex->second) = true;

          fv++;
        }
      }
      else
      {
        unsigned int fv = 0;
        for(auto iter = new_face_list.begin();
            iter != new_face_list.end();
            iter++)
        {
          new_face.vertices.at(fv) = iter->at(0);
          new_face.edges.at(fv == 0 ? (new_face_list.size() - 1) : (fv - 1)) = iter->at(1);

          const auto signed_vertex = graph.signed_vertices.find(iter->at(0));
          if (signed_vertex != graph.signed_vertices.end())
            vertices_included_in_face.at(signed_vertex->second) = true;

          fv++;
        }
      }
    }

    for (const auto& s_v : graph.signed_vertices)
    {
      if (!vertices_included_in_face.at(s_v.second))
        throw std::runtime_error("New faces not found for vertex " + std::to_string(s_v.first));
    }

    result.shrink_to_fit();

    return result;
  }
  // ***************************************************************************
  struct new_volume_data final
  {
      std::vector<unsigned int> vertices;
      std::vector<unsigned int> edges;
      std::vector<unsigned int> faces;
  };
  inline std::vector<new_volume_data> create_volumes_from_graph(const signed_graphs::graph& graph,
                                                                const std::vector<Gedim::MeshUtilities::MeshSurfaceIntersectionGraphData::FaceData>& graph_faces)
  {
    if (graph.vertices.size() < 4 ||
        graph.edges.size() < 6 ||
        graph.faces.size() < 4)
      return {};

    std::vector<new_volume_data> result;

    result.reserve(graph.vertices.size());

    std::vector<bool> faces_included_in_volume(graph.signed_faces.size(),
                                               false);

    for (const auto& face : graph.signed_faces)
    {
      if (faces_included_in_volume.at(face.second))
        continue;

      std::vector<bool> faces_visited(graph.faces.size(),
                                      false);

      for (const auto& s_f : graph.signed_faces)
      {
        const auto s_f_idx = graph.faces.at(s_f.first);
        faces_visited.at(s_f_idx) = faces_included_in_volume.at(s_f.second);
      }

      std::queue<unsigned int> queue;
      queue.push(face.first);
      while (!queue.empty())
      {
        const auto visited = queue.front();
        const auto visited_idx = graph.faces.at(visited);
        faces_visited.at(visited_idx) = true;
        queue.pop();

        for (const auto adjacent : graph.faces_edges_adjacency[visited_idx])
        {
          if (faces_visited.at(graph.faces.at(adjacent.face)))
            continue;

          queue.push(static_cast<int>(adjacent.face));
        }
      }

      std::vector<unsigned int> volume_faces;
      std::set<unsigned int> volume_vertices;
      std::set<unsigned int> volume_edges;

      volume_faces.reserve(faces_visited.size());
      for (const auto& v_f : graph.faces)
      {
        if (!faces_visited.at(v_f.second))
          continue;

        const auto signed_face = graph.signed_faces.find(v_f.first);
        if (signed_face != graph.signed_faces.end())
          faces_included_in_volume.at(signed_face->second) = true;

        volume_faces.push_back(v_f.first);
        const auto& graph_face = graph_faces.at(v_f.first);

        for (const auto f_v : graph_face.intersection_graph_vertices)
          volume_vertices.insert(f_v);
        for (const auto f_e : graph_face.intersection_graph_edges)
          volume_edges.insert(f_e);
      }
      volume_faces.shrink_to_fit();

      result.push_back({});
      auto& new_volume = result.back();
      new_volume.faces = volume_faces;
      new_volume.vertices.reserve(volume_vertices.size());
      new_volume.edges.reserve(volume_edges.size());

      for (const auto v_v : volume_vertices)
        new_volume.vertices.push_back(v_v);

      for (const auto v_e : volume_edges)
        new_volume.edges.push_back(v_e);
    }

    for (const auto& s_f : graph.signed_faces)
    {
      if (!faces_included_in_volume.at(s_f.second))
        throw std::runtime_error("New volumes not found for face " + std::to_string(s_f.first));
    }

    result.shrink_to_fit();

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

      std::map<unsigned int, unsigned int> original_face_vertex_to_mesh_cell1D;
      for (unsigned int fv = 0; fv < cell2D_vertices.size(); ++fv)
      {
        original_face_vertex_to_mesh_cell1D.insert(std::make_pair(result.mesh_cell0D_to_vertex.at(cell2D_vertices.at(fv)),
                                                                  cell2D_edges.at(fv)));
      }

      const auto& face_segments = surface_segments.cell2Ds_surface_segments.at(f);

      std::vector<unsigned int> face_vertices;
      std::map<unsigned int, unsigned int> face_edges;
      face_vertices.reserve(cell2D_vertices.size() +
                            cell2D_edges.size());

      for (const auto& cell0D_index : cell2D_vertices)
        face_vertices.push_back(result.mesh_cell0D_to_vertex.at(cell0D_index));

      for (const auto& cell1D_index : cell2D_edges)
      {
        const auto& surface_point_index = surface_intersections.cell1Ds_surface_point.at(cell1D_index);

        if (surface_point_index.index !=
            Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::None)
          face_vertices.push_back(result.surface_point_to_vertex.at(surface_point_index.index));

        for (const auto& edge_index : result.mesh_cell1D_to_edges.at(cell1D_index))
          face_edges.insert(std::make_pair(edge_index, face_edges.size()));
      }

      for (const auto& s : face_segments)
        face_edges.insert(std::make_pair(result.surface_segment_to_edge.at(s), face_edges.size()));

      face_vertices.shrink_to_fit();

      const auto face_graph = create_signed_graphs(face_vertices,
                                                   face_edges,
                                                   {},
                                                   result.vertices,
                                                   result.edges,
                                                   {});

      if (face_graph.is_only_negative ||
          face_graph.is_only_positive ||
          face_graph.is_only_on_surface)
      {
        const auto face_sign = face_graph.is_only_on_surface ?
                                 Gedim::MeshUtilities::PointSurfacePositions::OnSurface :
                                 face_graph.is_only_negative ?
                                   Gedim::MeshUtilities::PointSurfacePositions::Negative :
                                   Gedim::MeshUtilities::PointSurfacePositions::Positive;

        std::vector<unsigned int> face_edges_array;
        face_edges_array.reserve(face_edges.size());
        for (const auto& f_e : face_edges)
          face_edges_array.push_back(f_e.first);

        result.mesh_cell2D_to_faces.at(f) = { static_cast<unsigned int>(result.faces.size()) };
        result.faces.push_back(
              {
                face_vertices,
                face_edges_array,
                face_sign,
                f,
                Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::None
              });
        continue;
      }

      const auto new_negative_faces = create_faces_from_graph(original_face_vertex_to_mesh_cell1D,
                                                              face_graph.negative_graph,
                                                              result.edges);

      const auto new_positive_faces = create_faces_from_graph(original_face_vertex_to_mesh_cell1D,
                                                              face_graph.positive_graph,
                                                              result.edges);

      result.mesh_cell2D_to_faces.at(f).reserve(new_negative_faces.size() +
                                                new_positive_faces.size());

      for (const auto& new_face : new_negative_faces)
      {
        result.mesh_cell2D_to_faces.at(f).push_back(result.faces.size());
        result.faces.push_back(
              {
                new_face.vertices,
                new_face.edges,
                Gedim::MeshUtilities::PointSurfacePositions::Negative,
                f,
                Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::None
              });
      }

      for (const auto& new_face : new_positive_faces)
      {
        result.mesh_cell2D_to_faces.at(f).push_back(result.faces.size());
        result.faces.push_back(
              {
                new_face.vertices,
                new_face.edges,
                Gedim::MeshUtilities::PointSurfacePositions::Positive,
                f,
                Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::None
              });
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

    result.faces.shrink_to_fit();

    result.mesh_cell3D_to_volumes.resize(mesh.Cell3DTotalNumber());
    result.volumes.reserve(5 * mesh.Cell3DTotalNumber());

    for (unsigned int p = 0; p < mesh.Cell3DTotalNumber(); ++p)
    {
      const auto cell3D_vertices = mesh.Cell3DVertices(p);
      const auto cell3D_edges = mesh.Cell3DEdges(p);
      const auto cell3D_faces = mesh.Cell3DFaces(p);

      const auto& volume_triangles = surface_triangles.cell3Ds_surface_triangles.at(p);

      std::vector<unsigned int> volume_vertices;
      std::map<unsigned int, unsigned int> volume_edges;
      std::map<unsigned int, unsigned int> volume_faces;
      volume_vertices.reserve(cell3D_vertices.size() +
                              cell3D_edges.size());

      for (const auto& cell0D_index : cell3D_vertices)
        volume_vertices.push_back(result.mesh_cell0D_to_vertex.at(cell0D_index));

      for (const auto& cell1D_index : cell3D_edges)
      {
        const auto& surface_point_index = surface_intersections.cell1Ds_surface_point.at(cell1D_index);

        if (surface_point_index.index !=
            Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::None)
          volume_vertices.push_back(result.surface_point_to_vertex.at(surface_point_index.index));

        for (const auto& edge_index : result.mesh_cell1D_to_edges.at(cell1D_index))
          volume_edges.insert(std::make_pair(edge_index, volume_edges.size()));
      }

      for (const auto& cell2D_index : cell3D_faces)
      {
        for (const auto& s : surface_segments.cell2Ds_surface_segments.at(cell2D_index))
          volume_edges.insert(std::make_pair(result.surface_segment_to_edge.at(s), volume_edges.size()));

        for (const auto& face_index : result.mesh_cell2D_to_faces.at(cell2D_index))
          volume_faces.insert(std::make_pair(face_index, volume_faces.size()));
      }

      for (const auto& t : volume_triangles)
      {
        const auto& triangle_data = surface_triangles.surface_triangles.at(t);

        const std::vector<unsigned int> face_vertices = {
          result.surface_point_to_vertex.at(triangle_data.indices.at(0)),
          result.surface_point_to_vertex.at(triangle_data.indices.at(1)),
          result.surface_point_to_vertex.at(triangle_data.indices.at(2))
        };
        const auto face_edges = find_face_edges(face_vertices,
                                                result.edges_map);

        for (const auto t_e : face_edges)
        {
          if (!volume_edges.contains(t_e))
            volume_edges.insert(std::make_pair(t_e, volume_edges.size()));
        }

        volume_faces.insert(std::make_pair(result.surface_triangle_to_face.at(t), volume_faces.size()));
      }

      volume_vertices.shrink_to_fit();

      const auto volume_graph = create_signed_graphs(volume_vertices,
                                                     volume_edges,
                                                     volume_faces,
                                                     result.vertices,
                                                     result.edges,
                                                     result.faces);

      if (volume_graph.is_only_negative ||
          volume_graph.is_only_positive ||
          volume_graph.is_only_on_surface)
      {
        const auto volume_sign = volume_graph.is_only_on_surface ?
                                   Gedim::MeshUtilities::PointSurfacePositions::OnSurface :
                                   volume_graph.is_only_negative ?
                                     Gedim::MeshUtilities::PointSurfacePositions::Negative :
                                     Gedim::MeshUtilities::PointSurfacePositions::Positive;

        std::vector<unsigned int> volume_edges_array;
        volume_edges_array.reserve(volume_edges.size());
        for (const auto& v_e : volume_edges)
          volume_edges_array.push_back(v_e.first);

        std::vector<unsigned int> volume_faces_array;
        volume_faces_array.reserve(volume_faces.size());
        for (const auto& v_f : volume_faces)
          volume_faces_array.push_back(v_f.first);

        result.mesh_cell3D_to_volumes.at(p) = { static_cast<unsigned int>(result.volumes.size()) };
        result.volumes.push_back(
              {
                volume_vertices,
                volume_edges_array,
                volume_faces_array,
                volume_sign,
                p
              });

        continue;
      }

      const auto new_negative_volumes = create_volumes_from_graph(volume_graph.negative_graph,
                                                                  result.faces);

      const auto new_positive_volumes = create_volumes_from_graph(volume_graph.positive_graph,
                                                                  result.faces);

      result.mesh_cell3D_to_volumes.at(p).reserve(new_negative_volumes.size() +
                                                  new_positive_volumes.size());

      for (const auto& new_volume : new_negative_volumes)
      {
        result.mesh_cell3D_to_volumes.at(p).push_back(result.volumes.size());
        result.volumes.push_back(
              {
                new_volume.vertices,
                new_volume.edges,
                new_volume.faces,
                Gedim::MeshUtilities::PointSurfacePositions::Negative,
                p
              });
      }

      for (const auto& new_volume : new_positive_volumes)
      {
        result.mesh_cell3D_to_volumes.at(p).push_back(result.volumes.size());
        result.volumes.push_back(
              {
                new_volume.vertices,
                new_volume.edges,
                new_volume.faces,
                Gedim::MeshUtilities::PointSurfacePositions::Positive,
                p
              });
      }
    }

    result.volumes.shrink_to_fit();

    return result;
  }
  // ***************************************************************************
} // namespace Gedim
