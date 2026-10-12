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
Gedim::MeshUtilities::MarchingCubesData Gedim::MeshUtilities::MarchingCubes(
    const Gedim::GeometryUtilities &geometry_utilities,
    const std::function<Eigen::VectorXd(const Eigen::MatrixXd &)> &surface_implicit_equation,
    const Gedim::IMeshDAO &mesh)
{
    Gedim::MeshUtilities::MarchingCubesData result;

    auto marking_function = [&geometry_utilities,
                            &surface_implicit_equation](const Eigen::MatrixXd &points) {
        Eigen::VectorXi marks(points.cols());
        const auto surface_values = surface_implicit_equation(points);

        for (unsigned int p = 0; p < surface_values.size(); p++)
        {
            const auto surface_value = surface_values[p];

            if (geometry_utilities.IsValuePositive(surface_value,
                                                   geometry_utilities.Tolerance1D()))
            {
                marks[p] = +1.0;
            }
            else if (geometry_utilities.IsValueZero(surface_value,
                                                    geometry_utilities.Tolerance1D()))
              marks[p] = +0.0;
            else
              marks[p] = -1.0;
        }

        return marks;
    };

    std::vector<Eigen::MatrixXd> cell0Ds_coordinates;
    {
      const auto cell0Ds_coordinates_mat = mesh.Cell0DsCoordinates();

      cell0Ds_coordinates.resize(cell0Ds_coordinates_mat.cols());
      for (unsigned int c = 0; c < cell0Ds_coordinates_mat.cols(); ++c)
        cell0Ds_coordinates.at(c) = Eigen::Vector3d(cell0Ds_coordinates_mat.col(c));
    }

    result.cell0Ds_marked = MarkCells(marking_function,
                                      cell0Ds_coordinates,
                                      0.0);

    return result;
}
// ***************************************************************************
} // namespace Gedim
