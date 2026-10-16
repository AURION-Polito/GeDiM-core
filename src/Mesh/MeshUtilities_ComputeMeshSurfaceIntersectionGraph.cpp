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
  Gedim::MeshUtilities::MeshSurfaceIntersectionGraphData Gedim::MeshUtilities::ComputeMeshSurfaceIntersectionGraph(const std::vector<Gedim::MeshUtilities::PointSurfacePositions>& cell0Ds_surface_position,
                                                                                                                   const Gedim::MeshUtilities::IntersectMeshWithSurfaceResult& surface_intersections,
                                                                                                                   const Gedim::MeshUtilities::MarchingSquaresResult::SurfaceSegmentsData& surface_segments,
                                                                                                                   const Gedim::IMeshDAO &mesh)
  {
    Gedim::MeshUtilities::MeshSurfaceIntersectionGraphData result;
    
    result.vertices.reserve(mesh.Cell0DTotalNumber() +
                            surface_intersections.surface_intersection_points.size());

    for (unsigned int v = 0; v < mesh.Cell0DTotalNumber(); ++v)
    {
      result.vertices.push_back({
                                  mesh.Cell0DCoordinates(v),
                                  cell0Ds_surface_position.at(v),
                                  Gedim::MeshUtilities::MeshSurfaceIntersectionGraphData::VertexData::MeshTypes::Vertex,
                                  v,
                                  surface_intersections.cell0Ds_surface_points.at(v)
                                });
    }

    for (unsigned int v = 0; v < surface_intersections.surface_intersection_points.size(); ++v)
    {
      const auto& surface_intersection = surface_intersections.surface_intersection_points.at(v);

      switch (surface_intersection.Type)
      {
        case Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::SurfacePointData::Types::Vertex:
          continue;
        case Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::SurfacePointData::Types::Edge:
          break;
        default:
          throw std::runtime_error("Unknown intersection type");
      }

      result.vertices.push_back({
                                  surface_intersection.coordinate,
                                  Gedim::MeshUtilities::PointSurfacePositions::OnSurface,
                                  Gedim::MeshUtilities::MeshSurfaceIntersectionGraphData::VertexData::MeshTypes::Edge,
                                  surface_intersection.mesh_id,
                                  v
                                });
    }

    result.vertices.shrink_to_fit();

    return result;
  }
  // ***************************************************************************
} // namespace Gedim
