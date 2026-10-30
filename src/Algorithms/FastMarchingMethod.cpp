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

#include "FastMarchingMethod.hpp"
#include <queue>

namespace Gedim
{
  namespace Algorithms
  {
    namespace FastMarchingMethod
    {
      // ***************************************************************************
      Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodResult fast_marching_method_2D(const GeometryUtilities& geometry_utilities,
                                                                                              const Eigen::MatrixXd& mesh_vertices,
                                                                                              const std::vector<std::list<unsigned int> >& mesh_adj_list,
                                                                                              const double& constant_speed,
                                                                                              const unsigned int source_index)
      {
        using pq_item = std::pair<double, unsigned int>;
        using status_types = Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodResult::StatusTypes;

        Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodResult result;

        constexpr static double infty = std::numeric_limits<double>::infinity();

        result.vertices_status.resize(mesh_vertices.cols(), status_types::Far);
        result.vertices_time.resize(mesh_vertices.cols(), infty);

        result.vertices_status.at(source_index) = status_types::Verified;
        result.vertices_time.at(source_index) = 0.0;

        std::priority_queue<pq_item, std::vector<pq_item>, std::greater<pq_item>> min_pq;
        for (const auto& vn : mesh_adj_list.at(source_index))
        {
          result.vertices_time.at(vn) = geometry_utilities.SegmentLength(mesh_vertices.col(source_index),
                                                                         mesh_vertices.col(vn));
          result.vertices_status.at(source_index) = status_types::Trial;
          min_pq.push({ result.vertices_time.at(vn), vn });
        }

        return result;
      }
      // ***************************************************************************
      Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodResult fast_marching_method_2D(const GeometryUtilities& geometry_utilities,
                                                                                              const IMeshDAO& mesh,
                                                                                              const double& constant_speed,
                                                                                              const unsigned int source_index)
      {
        std::vector<std::list<unsigned int>> mesh_adj_list(mesh.Cell0DTotalNumber());
        for (unsigned int v = 0; v < mesh.Cell0DTotalNumber(); ++v)
        {
          if (!mesh.Cell0DIsActive(v))
            continue;

          for (unsigned int vn = 0; vn < mesh.Cell0DNumberNeighbourCell1D(v); ++vn)
          {
            if (!mesh.Cell0DHasNeighbourCell1D(v, vn))
              continue;

            mesh_adj_list.at(v).push_back(mesh.Cell0DNeighbourCell1D(v, vn));
          }
        }

        return fast_marching_method_2D(geometry_utilities,
                                       mesh.Cell0DsCoordinates(),
                                       mesh_adj_list,
                                       constant_speed,
                                       source_index);
      }
      // ***************************************************************************
    } // namespace FastMarchingMethod
  } // namespace Algorithms
} // namespace Gedim
