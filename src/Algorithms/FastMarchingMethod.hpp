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

#ifndef __Fast_Marching_Method_H
#define __Fast_Marching_Method_H

#include "Eigen/Eigen"
#include "GeometryUtilities.hpp"
#include "IMeshDAO.hpp"
#include <list>
#include <vector>

namespace Gedim
{
  namespace Algorithms
  {
    namespace FastMarchingMethod
    {
      struct FastMarchingMethodResult
      {
          enum struct StatusTypes
          {
            Far = 0,
            Trial = 1,
            Verified = 2
          };

          std::vector<double> vertices_time;
          std::vector<Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodResult::StatusTypes> vertices_status;
      };

      Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodResult fast_marching_method_2D(const Gedim::GeometryUtilities& geometry_utilities,
                                                                                              const Eigen::MatrixXd& mesh_vertices,
                                                                                              const std::vector<std::list<unsigned int>>& mesh_adj_list,
                                                                                              const double& constant_speed,
                                                                                              const unsigned int source_index);
      Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodResult fast_marching_method_2D(const Gedim::GeometryUtilities& geometry_utilities,
                                                                                              const Gedim::IMeshDAO& mesh,
                                                                                              const double& constant_speed,
                                                                                              const unsigned int source_index);

    } // namespace FastMarchingMethod
  } // namespace Algorithms
} // namespace Gedim

#endif // __Fast_Marching_Method_H
