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
    const std::function<Eigen::VectorXd(const Eigen::MatrixXd &)> &surface_implicit_equation,
    const Gedim::IMeshDAO &mesh)
{
    Gedim::Utilities::Unused(geometry_utilities);

    using SurfacePosition = Gedim::MeshUtilities::SurfacePositions;

    Gedim::MeshUtilities::IntersectMeshWithSurfaceResult result;

    result.cell0Ds_surface_point.resize(mesh.Cell0DTotalNumber());
    result.cell1Ds_surface_point.resize(mesh.Cell1DTotalNumber());

    result.surface_intersection_points.reserve(mesh.Cell0DTotalNumber() + mesh.Cell1DTotalNumber());

    const auto cell0Ds_surface_value = surface_implicit_equation(mesh.Cell0DsCoordinates());
    for (unsigned int p = 0; p < mesh.Cell0DTotalNumber(); p++)
    {
        auto &cell0D_surface_data = result.cell0Ds_surface_point.at(p);
        cell0D_surface_data.surface_value = cell0Ds_surface_value[p];
        cell0D_surface_data.index = Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::None;

        if (geometry_utilities.IsValuePositive(cell0D_surface_data.surface_value, geometry_utilities.Tolerance1D()))
            cell0D_surface_data.surface_position = SurfacePosition::Positive;
        else if (geometry_utilities.IsValueZero(cell0D_surface_data.surface_value, geometry_utilities.Tolerance1D()))
            cell0D_surface_data.surface_position = SurfacePosition::OnSurface;
        else if (geometry_utilities.IsValueNegative(cell0D_surface_data.surface_value, geometry_utilities.Tolerance1D()))
            cell0D_surface_data.surface_position = SurfacePosition::Negative;
        else
            throw std::runtime_error("unknown cell0D position");

        switch (cell0D_surface_data.surface_position)
        {
        case SurfacePosition::OnSurface: {
            cell0D_surface_data.index = result.surface_intersection_points.size();
            result.surface_intersection_points.push_back(
                {mesh.Cell0DCoordinates(p), Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::SurfacePointData::Types::Vertex, p, 0.0});
        }
        break;
        case SurfacePosition::Positive:
        case SurfacePosition::Negative:
            break;
        default:
            throw std::runtime_error("Position not managed");
        }
    }

    for (unsigned int e = 0; e < mesh.Cell1DTotalNumber(); ++e)
    {
        auto &cell1D_surface_data = result.cell1Ds_surface_point.at(e);
        cell1D_surface_data.surface_value = 0.0;
        cell1D_surface_data.index = Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::None;

        const auto origin_id = mesh.Cell1DOrigin(e);
        const auto end_id = mesh.Cell1DEnd(e);

        const auto origin_coordinate = mesh.Cell1DOriginCoordinates(e);
        const auto end_coordinate = mesh.Cell1DEndCoordinates(e);

        const auto origin_surface_position = result.cell0Ds_surface_point.at(origin_id);
        const auto end_surface_position = result.cell0Ds_surface_point.at(end_id);

        switch (origin_surface_position.surface_position)
        {
        case SurfacePosition::OnSurface: {
            switch (end_surface_position.surface_position)
            {
            case SurfacePosition::OnSurface: {
                cell1D_surface_data.surface_position = SurfacePosition::OnSurface;
            }
                continue;
            case SurfacePosition::Positive:
            case SurfacePosition::Negative: {
                cell1D_surface_data.surface_position = end_surface_position.surface_position;
            }
                continue;
            default:
                throw std::runtime_error("Position not managed");
            }
        }
        break;
        case SurfacePosition::Positive:
        case SurfacePosition::Negative: {
            switch (end_surface_position.surface_position)
            {
            case SurfacePosition::OnSurface: {
                cell1D_surface_data.surface_position = origin_surface_position.surface_position;
            }
                continue;
            case SurfacePosition::Positive:
            case SurfacePosition::Negative:
                break;
            default:
                throw std::runtime_error("Position not managed");
            }
        }
        break;
        default:
            throw std::runtime_error("Position not managed");
        }

        if (origin_surface_position.surface_position != end_surface_position.surface_position)
        {
            const double origin_surface_value = cell0Ds_surface_value[origin_id];
            const double end_surface_value = cell0Ds_surface_value[end_id];

            const double curvilinear_coordinate = origin_surface_value / (origin_surface_value - end_surface_value);

            const Eigen::Vector3d intersection = origin_coordinate + curvilinear_coordinate * (end_coordinate - origin_coordinate);
            const auto cell1Ds_surface_value = surface_implicit_equation(intersection);

            cell1D_surface_data.surface_value = cell1Ds_surface_value[0];
            cell1D_surface_data.index = result.surface_intersection_points.size();
            cell1D_surface_data.surface_position = SurfacePosition::OnSurface;
            result.surface_intersection_points.push_back(
                {intersection, Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::SurfacePointData::Types::Edge, e, curvilinear_coordinate});
        }
        else
            cell1D_surface_data.surface_position = origin_surface_position.surface_position;
    }

    result.surface_intersection_points.shrink_to_fit();

    return result;
}
// ***************************************************************************
} // namespace Gedim
