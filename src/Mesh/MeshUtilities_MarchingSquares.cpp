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

#include "CommonUtilities.hpp"
#include "MeshUtilities.hpp"

namespace Gedim
{
  // ***************************************************************************
  inline MeshUtilities::MarchingSquaresResult::SurfaceTrianglesData CreateTriangles(
      const Gedim::GeometryUtilities &geometry_utilities,
      const Eigen::VectorXd &cell0Ds_surface_value,
      const std::vector<Gedim::MeshUtilities::PointSurfacePositions> &cell0Ds_surface_position,
      const std::vector<Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::SurfacePointData> &surface_intersection_points,
      const std::vector<unsigned int> &cell0Ds_surface_points,
      const std::vector<unsigned int> &cell1Ds_surface_points,
      const Gedim::IMeshDAO &mesh)
  {
    MeshUtilities::MarchingSquaresResult::SurfaceTrianglesData result;

    return result;
  }
  // ***************************************************************************
  Gedim::MeshUtilities::MarchingSquaresResult Gedim::MeshUtilities::MarchingSquares(
      const Gedim::GeometryUtilities &geometry_utilities,
      const std::function<Eigen::VectorXd(const Eigen::MatrixXd &)> &surface_implicit_equation,
      const Gedim::IMeshDAO &mesh)
  {
    using SurfacePosition = Gedim::MeshUtilities::PointSurfacePositions;

    Gedim::MeshUtilities::MarchingSquaresResult result;

    result.cell0Ds_surface_value = surface_implicit_equation(mesh.Cell0DsCoordinates());

    {
      auto marking_function = [&geometry_utilities](const Eigen::VectorXd &points_surface_value) {
        std::vector<SurfacePosition> marks(points_surface_value.size());

        for (unsigned int p = 0; p < points_surface_value.size(); p++)
        {
          const auto surface_value = points_surface_value[p];

          if (geometry_utilities.IsValuePositive(surface_value, geometry_utilities.Tolerance1D()))
          {
            marks[p] = SurfacePosition::Positive;
          }
          else if (geometry_utilities.IsValueZero(surface_value, geometry_utilities.Tolerance1D()))
            marks[p] = SurfacePosition::OnSurface;
          else if (geometry_utilities.IsValueNegative(surface_value, geometry_utilities.Tolerance1D()))
            marks[p] = SurfacePosition::Negative;
          else
            throw std::runtime_error("unknown cell0D position");
        }

        return marks;
      };

      result.cell0Ds_surface_position = marking_function(result.cell0Ds_surface_value);
    }

    result.surface_intersections =
        IntersectMeshWithSurface(geometry_utilities, result.cell0Ds_surface_value, result.cell0Ds_surface_position, mesh);

    result.surface_triangles = CreateTriangles(geometry_utilities,
                                               result.cell0Ds_surface_value,
                                               result.cell0Ds_surface_position,
                                               result.surface_intersections.surface_intersection_points,
                                               result.surface_intersections.cell0Ds_surface_points,
                                               result.surface_intersections.cell1Ds_surface_points,
                                               mesh);

    return result;
  }
  // ***************************************************************************
} // namespace Gedim
