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
  Gedim::MeshUtilities::IntersectMeshWithSurfaceResult Gedim::MeshUtilities::IntersectMeshWithSurface(const Gedim::GeometryUtilities &geometry_utilities,
                                                                                                      const Eigen::VectorXd &cell0Ds_surface_value,
                                                                                                      const std::vector<Gedim::MeshUtilities::MarchingCubesResult::SurfacePositions> &cell0Ds_surface_position,
                                                                                                      const Gedim::IMeshDAO &mesh)
  {
    Gedim::MeshUtilities::IntersectMeshWithSurfaceResult result;

    result.cell0Ds_surface_points.resize(mesh.Cell0DTotalNumber(),
                                         Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::None);
    result.cell1Ds_surface_points.resize(mesh.Cell1DTotalNumber(),
                                         Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::None);
    result.surface_intersection_points.reserve(mesh.Cell1DTotalNumber() +
                                               mesh.Cell0DTotalNumber());

    for (unsigned int e = 0; e < mesh.Cell1DTotalNumber(); ++e)
    {
      const auto origin_id = mesh.Cell1DOrigin(e);
      const auto end_id = mesh.Cell1DEnd(e);

      const auto origin_coordinate = mesh.Cell1DOriginCoordinates(e);
      const auto end_coordinate = mesh.Cell1DEndCoordinates(e);

      const auto origin_surface_position = cell0Ds_surface_position.at(origin_id);
      const auto end_surface_position = cell0Ds_surface_position.at(end_id);

      switch (origin_surface_position)
      {
        case Gedim::MeshUtilities::MarchingCubesResult::SurfacePositions::OnSurface:
        {
          if (result.cell0Ds_surface_points.at(origin_id) ==
              Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::None)
          {
            result.cell0Ds_surface_points.at(origin_id) =
                result.surface_intersection_points.size();
            result.surface_intersection_points.push_back({
                                                           origin_coordinate,
                                                           Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::SurfacePointData::Types::Vertex,
                                                           origin_id
                                                         });
          }
        }
          break;
        case Gedim::MeshUtilities::MarchingCubesResult::SurfacePositions::Positive:
        case Gedim::MeshUtilities::MarchingCubesResult::SurfacePositions::Negative:
          break;
        default:
          throw std::runtime_error("Position not managed");
      }

      switch (end_surface_position)
      {
        case Gedim::MeshUtilities::MarchingCubesResult::SurfacePositions::OnSurface:
        {
          if (result.cell0Ds_surface_points.at(end_id) ==
              Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::None)
          {
            result.cell0Ds_surface_points.at(end_id) =
                result.surface_intersection_points.size();
            result.surface_intersection_points.push_back({
                                                           end_coordinate,
                                                           Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::SurfacePointData::Types::Vertex,
                                                           end_id
                                                         });
          }
        }
          break;
        case Gedim::MeshUtilities::MarchingCubesResult::SurfacePositions::Positive:
        case Gedim::MeshUtilities::MarchingCubesResult::SurfacePositions::Negative:
          break;
        default:
          throw std::runtime_error("Position not managed");
      }

      if (origin_surface_position != end_surface_position &&
          origin_surface_position != Gedim::MeshUtilities::MarchingCubesResult::SurfacePositions::OnSurface &&
          end_surface_position != Gedim::MeshUtilities::MarchingCubesResult::SurfacePositions::OnSurface)
      {
        const double origin_surface_value = cell0Ds_surface_value[origin_id];
        const double end_surface_value = cell0Ds_surface_value[end_id];

        const Eigen::Vector3d intersection = origin_coordinate +
                                             origin_surface_value *
                                             (end_coordinate - origin_coordinate) /
                                             (origin_surface_value - end_surface_value);

        result.cell1Ds_surface_points.at(e) =
            result.surface_intersection_points.size();
        result.surface_intersection_points.push_back({
                                                       intersection,
                                                       Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::SurfacePointData::Types::Edge,
                                                       e
                                                     });

      }
    }

    return result;
  }
  // ***************************************************************************
  Gedim::MeshUtilities::MarchingCubesResult Gedim::MeshUtilities::MarchingCubes(
      const Gedim::GeometryUtilities &geometry_utilities,
      const std::function<Eigen::VectorXd(const Eigen::MatrixXd &)> &surface_implicit_equation,
      const Gedim::IMeshDAO &mesh)
  {
    Gedim::MeshUtilities::MarchingCubesResult result;

    result.cell0Ds_surface_value = surface_implicit_equation(mesh.Cell0DsCoordinates());

    {
      auto marking_function = [&geometry_utilities](const Eigen::VectorXd &points_surface_value) {
        std::vector<Gedim::MeshUtilities::MarchingCubesResult::SurfacePositions> marks(points_surface_value.size());

        for (unsigned int p = 0; p < points_surface_value.size(); p++)
        {
          const auto surface_value = points_surface_value[p];

          if (geometry_utilities.IsValuePositive(surface_value, geometry_utilities.Tolerance1D()))
          {
            marks[p] = Gedim::MeshUtilities::MarchingCubesResult::SurfacePositions::Positive;
          }
          else if (geometry_utilities.IsValueZero(surface_value, geometry_utilities.Tolerance1D()))
            marks[p] = Gedim::MeshUtilities::MarchingCubesResult::SurfacePositions::OnSurface;
          else if (geometry_utilities.IsValueNegative(surface_value, geometry_utilities.Tolerance1D()))
            marks[p] = Gedim::MeshUtilities::MarchingCubesResult::SurfacePositions::Negative;
          else
            throw std::runtime_error("unknown cell0D position");
        }

        return marks;
      };

      result.cell0Ds_surface_position = marking_function(result.cell0Ds_surface_value);
    }

    result.surface_intersections = IntersectMeshWithSurface(geometry_utilities,
                                                            result.cell0Ds_surface_value,
                                                            result.cell0Ds_surface_position,
                                                            mesh);

    return result;
  }
  // ***************************************************************************
} // namespace Gedim
