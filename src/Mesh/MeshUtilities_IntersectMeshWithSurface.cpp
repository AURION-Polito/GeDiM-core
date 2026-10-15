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
Gedim::MeshUtilities::IntersectMeshWithSurfaceResult Gedim::MeshUtilities::IntersectMeshWithSurface(
    const Gedim::GeometryUtilities &geometry_utilities,
    const Eigen::VectorXd &cell0Ds_surface_value,
    const std::vector<Gedim::MeshUtilities::PointSurfacePositions> &cell0Ds_surface_position,
    const Gedim::IMeshDAO &mesh)
{
    Gedim::Utilities::Unused(geometry_utilities);

    using SurfacePosition = Gedim::MeshUtilities::PointSurfacePositions;

    Gedim::MeshUtilities::IntersectMeshWithSurfaceResult result;

    result.cell0Ds_surface_points.resize(mesh.Cell0DTotalNumber(), Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::None);
    result.cell1Ds_surface_points.resize(mesh.Cell1DTotalNumber(), Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::None);
    result.surface_intersection_points.reserve(mesh.Cell1DTotalNumber() + mesh.Cell0DTotalNumber());

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
        case SurfacePosition::OnSurface: {
            if (result.cell0Ds_surface_points.at(origin_id) == Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::None)
            {
                result.cell0Ds_surface_points.at(origin_id) = result.surface_intersection_points.size();
                result.surface_intersection_points.push_back(
                    {origin_coordinate, Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::SurfacePointData::Types::Vertex, origin_id});
            }
        }
        break;
        case SurfacePosition::Positive:
        case SurfacePosition::Negative:
            break;
        default:
            throw std::runtime_error("Position not managed");
        }

        switch (end_surface_position)
        {
        case SurfacePosition::OnSurface: {
            if (result.cell0Ds_surface_points.at(end_id) == Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::None)
            {
                result.cell0Ds_surface_points.at(end_id) = result.surface_intersection_points.size();
                result.surface_intersection_points.push_back(
                    {end_coordinate, Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::SurfacePointData::Types::Vertex, end_id});
            }
        }
        break;
        case SurfacePosition::Positive:
        case SurfacePosition::Negative:
            break;
        default:
            throw std::runtime_error("Position not managed");
        }

        if (origin_surface_position != end_surface_position && origin_surface_position != SurfacePosition::OnSurface &&
            end_surface_position != SurfacePosition::OnSurface)
        {
            const double origin_surface_value = cell0Ds_surface_value[origin_id];
            const double end_surface_value = cell0Ds_surface_value[end_id];

            const Eigen::Vector3d intersection = origin_coordinate + origin_surface_value * (end_coordinate - origin_coordinate) /
                                                                         (origin_surface_value - end_surface_value);

            result.cell1Ds_surface_points.at(e) = result.surface_intersection_points.size();
            result.surface_intersection_points.push_back(
                {intersection, Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::SurfacePointData::Types::Edge, e});
        }
    }

    result.surface_intersection_points.shrink_to_fit();

    return result;
}
// ***************************************************************************
} // namespace Gedim
