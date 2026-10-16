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
  Gedim::MeshUtilities::MeshSurfaceIntersectionGraphData Gedim::MeshUtilities::ComputeMeshSurfaceIntersectionGraph(const Gedim::MeshUtilities::IntersectMeshWithSurfaceResult& surface_intersections,
                                                                                                                   const Gedim::MeshUtilities::MarchingSquaresResult::SurfaceSegmentsData& surface_segments,
                                                                                                                   const Gedim::IMeshDAO &mesh)
  {
    Gedim::MeshUtilities::MeshSurfaceIntersectionGraphData result;
    
    result.mesh_cell0D_to_vertex.resize(mesh.Cell0DTotalNumber());
    result.surface_point_to_vertex.resize(surface_intersections.surface_intersection_points.size());
    result.vertices.reserve(mesh.Cell0DTotalNumber() +
                            surface_intersections.surface_intersection_points.size());

    for (unsigned int v = 0; v < mesh.Cell0DTotalNumber(); ++v)
    {
      const auto vertex_surface_point = surface_intersections.cell0Ds_surface_point.at(v);

      switch (vertex_surface_point.surface_position)
      {
        case Gedim::MeshUtilities::PointSurfacePositions::Positive:
        case Gedim::MeshUtilities::PointSurfacePositions::Negative:
        {
          result.mesh_cell0D_to_vertex.at(v) = result.vertices.size();
          result.vertices.push_back({
                                      mesh.Cell0DCoordinates(v),
                                      surface_intersections.cell0Ds_surface_point.at(v).surface_position,
                                      Gedim::MeshUtilities::MeshSurfaceIntersectionGraphData::VertexData::MeshTypes::Vertex,
                                      v,
                                      surface_intersections.cell0Ds_surface_point.at(v).index
                                    });
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
      const auto& surface_intersection = surface_intersections.surface_intersection_points.at(s);

      switch (surface_intersection.Type)
      {
        case Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::SurfacePointData::Types::Vertex:
        {
          result.mesh_cell0D_to_vertex.at(surface_intersection.mesh_id) = result.vertices.size();
          result.surface_point_to_vertex.at(s) = result.vertices.size();
          result.vertices.push_back({
                                      surface_intersection.coordinate,
                                      Gedim::MeshUtilities::PointSurfacePositions::OnSurface,
                                      Gedim::MeshUtilities::MeshSurfaceIntersectionGraphData::VertexData::MeshTypes::Vertex,
                                      surface_intersection.mesh_id,
                                      s
                                    });
        }
          continue;
        case Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::SurfacePointData::Types::Edge:
        {
          result.surface_point_to_vertex.at(s) = result.vertices.size();
          result.vertices.push_back({
                                      surface_intersection.coordinate,
                                      Gedim::MeshUtilities::PointSurfacePositions::OnSurface,
                                      Gedim::MeshUtilities::MeshSurfaceIntersectionGraphData::VertexData::MeshTypes::Edge,
                                      surface_intersection.mesh_id,
                                      s
                                    });
        }
          break;
        default:
          throw std::runtime_error("Unknown intersection type");
      }
    }

    result.vertices.shrink_to_fit();

    result.mesh_cell1D_to_edges.resize(mesh.Cell1DTotalNumber());
    result.surface_segment_to_edge.resize(surface_segments.surface_segments.size());
    result.edges.reserve(2 * mesh.Cell1DTotalNumber() +
                         surface_segments.surface_segments.size());

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
          result.mesh_cell1D_to_edges.at(e) = { result.edges.size() };
          result.edges.push_back({
                                   { edge_origin_id, edge_end_id },
                                   edge_surface_point.surface_position,
                                   e,
                                   Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::None
                                 });
        }
          continue;
        case Gedim::MeshUtilities::PointSurfacePositions::OnSurface:
        {
          const auto edge_origin_surface_position = surface_intersections.cell0Ds_surface_point.at(edge_origin_id).surface_position;
          const auto edge_end_surface_position = surface_intersections.cell0Ds_surface_point.at(edge_end_id).surface_position;

          if (edge_origin_surface_position == Gedim::MeshUtilities::PointSurfacePositions::OnSurface &&
              edge_end_surface_position == Gedim::MeshUtilities::PointSurfacePositions::OnSurface)
            continue;


        }
          break;
        default:
          throw std::runtime_error("Unknown position type");
      }
    }

    for (unsigned int s = 0; s < surface_segments.surface_segments.size(); ++s)
    {
      const auto& segment_data = surface_segments.surface_segments.at(s);

      result.surface_segment_to_edge.at(s) = result.edges.size();

      result.edges.push_back({
                               {
                                 result.surface_point_to_vertex.at(segment_data.indices.at(0)),
                                 result.surface_point_to_vertex.at(segment_data.indices.at(1))
                               },
                               Gedim::MeshUtilities::PointSurfacePositions::OnSurface,
                               Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::None,
                               s
                             });
    }

    result.edges.shrink_to_fit();

    return result;
  }
  // ***************************************************************************
} // namespace Gedim
