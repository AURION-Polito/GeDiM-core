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
      constexpr static double infty = std::numeric_limits<double>::infinity();
      // ***************************************************************************
      inline double compute_time_ij(const GeometryUtilities& geometry_utilities,
                                    const Eigen::MatrixXd& mesh_vertices,
                                    const double& constant_speed,
                                    const unsigned int i,
                                    const unsigned int j)
      {
        return geometry_utilities.IsValueZero(constant_speed,
                                              geometry_utilities.Tolerance1D()) ?
              infty :
              geometry_utilities.SegmentLength(mesh_vertices.col(i),
                                               mesh_vertices.col(j)) / constant_speed;
      }
      // ***************************************************************************
      inline double eikonal_triangle(const GeometryUtilities& geometry_utilities,
                                     const Eigen::MatrixXd& mesh_vertices,
                                     const double& constant_speed,
                                     const unsigned int i,
                                     const unsigned int j,
                                     const unsigned int k)
      {
        return std::min(compute_time_ij(geometry_utilities,
                                        mesh_vertices,
                                        constant_speed,
                                        i, k),
                        compute_time_ij(geometry_utilities,
                                        mesh_vertices,
                                        constant_speed,
                                        j, k));
      }
      // ***************************************************************************
      Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodResult fast_marching_method_2D(const GeometryUtilities& geometry_utilities,
                                                                                              const Eigen::MatrixXd& mesh_vertices,
                                                                                              const std::vector<std::list<unsigned int> >& mesh_adj_list,
                                                                                              const std::map<std::pair<unsigned int, unsigned int>, std::list<unsigned int>>& mesh_triangles,
                                                                                              const double& constant_speed,
                                                                                              const unsigned int source_index)
      {
        using pq_item = std::pair<double, unsigned int>;
        using status_types = Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodResult::StatusTypes;

        Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodResult result;

        result.vertices_status.resize(mesh_vertices.cols(), status_types::Far);
        result.vertices_time.resize(mesh_vertices.cols(), infty);

        result.vertices_status.at(source_index) = status_types::Verified;
        result.vertices_time.at(source_index) = 0.0;

        std::priority_queue<pq_item, std::vector<pq_item>, std::greater<pq_item>> min_pq;
        for (const auto& vn : mesh_adj_list.at(source_index))
        {
          result.vertices_time.at(vn) = result.vertices_time.at(source_index) +
                                        compute_time_ij(geometry_utilities,
                                                        mesh_vertices,
                                                        constant_speed,
                                                        source_index,
                                                        vn);

          result.vertices_status.at(vn) = status_types::Trial;
          min_pq.push({ result.vertices_time.at(vn), vn });
        }

        while (!min_pq.empty())
        {
          const auto selected_v = min_pq.top();
          min_pq.pop();

          if (result.vertices_status.at(selected_v.second) ==
              status_types::Verified)
            continue;

          result.vertices_status.at(selected_v.second) = status_types::Verified;

          for (const auto& adj_i : mesh_adj_list.at(selected_v.second))
          {
            if (result.vertices_status.at(adj_i) ==
                status_types::Verified)
              continue;

            auto new_time_i = std::min(result.vertices_time.at(adj_i),
                                       result.vertices_time.at(selected_v.second) +
                                       compute_time_ij(geometry_utilities,
                                                       mesh_vertices,
                                                       constant_speed,
                                                       selected_v.second,
                                                       adj_i));

            for (const auto k : mesh_triangles.at(std::make_pair(selected_v.second, adj_i)))
            {
              if (result.vertices_status.at(selected_v.second) !=
                  status_types::Verified)
                continue;

              const auto time_k = eikonal_triangle(geometry_utilities,
                                                   mesh_vertices,
                                                   constant_speed,
                                                   selected_v.second,
                                                   k,
                                                   adj_i);

              if (time_k < new_time_i)
                new_time_i = time_k;
            }

            if (new_time_i < result.vertices_time.at(adj_i))
            {
              result.vertices_time.at(adj_i) = new_time_i;
              result.vertices_status.at(adj_i) = status_types::Trial;
              min_pq.push({ result.vertices_time.at(adj_i), adj_i });
            }
          }
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
        std::map<std::pair<unsigned int, unsigned int>, std::list<unsigned int> > mesh_triangles;

        for (unsigned int v = 0; v < mesh.Cell0DTotalNumber(); ++v)
        {
          if (!mesh.Cell0DIsActive(v))
            continue;

          for (unsigned int vn = 0; vn < mesh.Cell0DNumberNeighbourCell1D(v); ++vn)
          {
            if (!mesh.Cell0DHasNeighbourCell1D(v, vn))
              continue;

            const auto cell1D = mesh.Cell0DNeighbourCell1D(v, vn);

            if (mesh.Cell1DOrigin(cell1D) == v)
              mesh_adj_list.at(v).push_back(mesh.Cell1DEnd(cell1D));
            else
              mesh_adj_list.at(v).push_back(mesh.Cell1DOrigin(cell1D));
          }
        }

        for (unsigned int f = 0; f < mesh.Cell2DTotalNumber(); ++f)
        {
          if (!mesh.Cell2DIsActive(f))
            continue;

          for (unsigned int fv = 0; fv < mesh.Cell2DNumberVertices(f); ++fv)
          {
            const unsigned int fv_i = mesh.Cell2DVertex(f, fv);
            const unsigned int fv_j = mesh.Cell2DVertex(f, (fv + 1) % mesh.Cell2DNumberVertices(f));
            const unsigned int fv_k = mesh.Cell2DVertex(f, (fv + 2) % mesh.Cell2DNumberVertices(f));

            if (!mesh_triangles.contains(std::make_pair(fv_i, fv_j)))
            {
              mesh_triangles.insert(std::make_pair(std::make_pair(fv_i, fv_j), std::list<unsigned int>()));
              mesh_triangles.insert(std::make_pair(std::make_pair(fv_j, fv_i), std::list<unsigned int>()));
            }

            mesh_triangles.at(std::make_pair(fv_i, fv_j)).push_back(fv_k);
            mesh_triangles.at(std::make_pair(fv_j, fv_i)).push_back(fv_k);
          }
        }

        return fast_marching_method_2D(geometry_utilities,
                                       mesh.Cell0DsCoordinates(),
                                       mesh_adj_list,
                                       mesh_triangles,
                                       constant_speed,
                                       source_index);
      }
      // ***************************************************************************
    } // namespace FastMarchingMethod
  } // namespace Algorithms
} // namespace Gedim
