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
#include "BrentSegmentMinimization.hpp"
#include <queue>

namespace Gedim
{
  namespace Algorithms
  {
    namespace FastMarchingMethod
    {
      // ***************************************************************************
      constexpr static double infty = std::numeric_limits<double>::infinity();
      using vertex_type = Gedim::Algorithms::FastMarchingMethod::MeshGraph::Vertex::VertexTypes;
      using vertex_initialization_types = Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::VertexInitialCondition::InitializationTypes;
      using edge_restriction_types = Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::EdgeRestriction::RestrictionTypes;
      using speed_types = Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::SpeedTypes;
      using eikonal_min_types = Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::EikonalMinimizationTypes;
      // ***************************************************************************
      inline double compute_time_ij(const GeometryUtilities& geometry_utilities,
                                    const double& constant_speed,
                                    const Eigen::Vector3d& x_i,
                                    const Eigen::Vector3d& x_j)
      {
        return geometry_utilities.IsValueZero(constant_speed,
                                              geometry_utilities.Tolerance1D()) ?
              infty :
              std::isinf(constant_speed) ?
                0.0 :
                geometry_utilities.SegmentLength(x_i,
                                                 x_j) / constant_speed;
      }
      // ***************************************************************************
      inline double eikonal_triangle(const GeometryUtilities& geometry_utilities,
                                     const MeshGraph& mesh_graph,
                                     const double& constant_speed,
                                     const std::vector<double>& vertices_time,
                                     const unsigned int i,
                                     const unsigned int j,
                                     const unsigned int k,
                                     const Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::EdgeRestriction& edge_ij_restriction,
                                     const Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::EikonalMinimizationTypes& miminization_type)
      {
        unsigned int e_ij_o = std::min(i, j);
        unsigned int e_ij_e = std::max(i, j);

        auto min_f =
            [&](double lambda)
        {
          return lambda * vertices_time.at(e_ij_o) +
              (1.0 - lambda) * vertices_time.at(e_ij_e) +
              compute_time_ij(geometry_utilities,
                              constant_speed,
                              mesh_graph.mesh_vertices.at(k).coordinate,
                              lambda * mesh_graph.mesh_vertices.at(e_ij_o).coordinate +
                              (1.0 - lambda) * mesh_graph.mesh_vertices.at(e_ij_e).coordinate);
        };

        std::array<double, 2> edge_reachable_zone;
        switch (edge_ij_restriction.type)
        {
          case edge_restriction_types::Unreachable:
            return infty;
          case edge_restriction_types::Free:
            edge_reachable_zone = { 0.0, 1.0 };
            break;
          case edge_restriction_types::Constrained:
          {
            edge_reachable_zone = edge_ij_restriction.reachable_zone;
          }
            break;
          default:
            throw std::runtime_error("unknown edge restriction type");
        }

        switch (miminization_type)
        {
          case eikonal_min_types::Simple:
          {
            return std::min({
                              min_f(edge_reachable_zone.at(0)),
                              min_f(0.5 * (edge_reachable_zone.at(1) + edge_reachable_zone.at(0))),
                              min_f(edge_reachable_zone.at(1))
                            });
          }
          case eikonal_min_types::Linear:
          {
            if (geometry_utilities.IsValueZero(std::abs(edge_reachable_zone.at(1) - edge_reachable_zone.at(0)),
                                               geometry_utilities.Tolerance1D()))
              return min_f(edge_reachable_zone.at(0));

            const auto min_result = Gedim::Algorithms::Minimization::brent_segment_minimization(
                                      min_f,
                                      edge_reachable_zone.at(0),
                                      edge_reachable_zone.at(1),
                                      1.0e-6,
                                      500
                                      );

            if (!min_result.converged)
            {
              //          std::cout.precision(2);
              //          std::cout<< std::scientific<< "V_i "<< min_f(0.0)<< " ";
              //          std::cout<< std::scientific<< "V_j "<< min_f(1.0)<< " ";
              //          std::cout<< std::scientific<< "k "<< min_result.min<< " ";
              //          std::cout<< std::scientific<< "V_k "<< min_result.f_min<< " ";

              //          Gedim::Algorithms::Minimization::brent_minimize(
              //                min_f,
              //                0.0,
              //                1.0,
              //                1.0e-6,
              //                500);

              throw std::runtime_error("Minimization algorithm failed");
            }

            return min_result.f_min;
          }
          default:
            throw std::runtime_error("Minimization type unknown");
        }
      }
      // ***************************************************************************
      MeshGraph create_mesh_graph(const IMeshDAO& mesh)
      {
        Gedim::Algorithms::FastMarchingMethod::MeshGraph mesh_graph;

        mesh_graph.mesh_vertices.resize(mesh.Cell0DTotalNumber());

        for (unsigned int v = 0; v < mesh.Cell0DTotalNumber(); ++v)
        {
          auto& mesh_vertex = mesh_graph.mesh_vertices.at(v);

          mesh_vertex.coordinate = mesh.Cell0DCoordinates(v);

          if (mesh.Cell0DMarker(v) == 0)
            mesh_vertex.type = vertex_type::Internal;
          else
            mesh_vertex.type = vertex_type::Border;

          for (unsigned int vn = 0; vn < mesh.Cell0DNumberNeighbourCell1D(v); ++vn)
          {
            if (!mesh.Cell0DHasNeighbourCell1D(v, vn))
              continue;

            const auto cell1D = mesh.Cell0DNeighbourCell1D(v, vn);

            if (mesh.Cell1DOrigin(cell1D) == v)
              mesh_vertex.adj_list.push_back(mesh.Cell1DEnd(cell1D));
            else
              mesh_vertex.adj_list.push_back(mesh.Cell1DOrigin(cell1D));
          }
        }

        for (unsigned int e = 0; e < mesh.Cell1DTotalNumber(); ++e)
        {

          const auto extremes = mesh.Cell1DExtremes(e);

          mesh_graph.mesh_edges.insert(std::make_pair(std::make_pair(extremes[0], extremes[1]), e));
          mesh_graph.mesh_edges.insert(std::make_pair(std::make_pair(extremes[1], extremes[0]), e));
        }

        for (unsigned int f = 0; f < mesh.Cell2DTotalNumber(); ++f)
        {
          for (unsigned int fv = 0; fv < mesh.Cell2DNumberVertices(f); ++fv)
          {
            const unsigned int fv_i = mesh.Cell2DVertex(f, fv);
            const unsigned int fv_j = mesh.Cell2DVertex(f, (fv + 1) % mesh.Cell2DNumberVertices(f));
            const unsigned int fv_k = mesh.Cell2DVertex(f, (fv + 2) % mesh.Cell2DNumberVertices(f));

            if (!mesh_graph.mesh_triangles.contains(std::make_pair(fv_i, fv_j)))
            {
              mesh_graph.mesh_triangles.insert(std::make_pair(std::make_pair(fv_i, fv_j), std::list<unsigned int>()));
              mesh_graph.mesh_triangles.insert(std::make_pair(std::make_pair(fv_j, fv_i), std::list<unsigned int>()));
            }

            mesh_graph.mesh_triangles.at(std::make_pair(fv_i, fv_j)).push_back(fv_k);
            mesh_graph.mesh_triangles.at(std::make_pair(fv_j, fv_i)).push_back(fv_k);
          }
        }

        return mesh_graph;
      }
      // ***************************************************************************
      Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig create_configuration(const GeometryUtilities& geometry_utilities,
                                                                                           const MeshGraph& mesh_graph,
                                                                                           const std::vector<unsigned int>& sources_index,
                                                                                           const double& ball_radius,
                                                                                           const double& constant_speed,
                                                                                           const Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::EikonalMinimizationTypes& eikonal_minimization_type)
      {
        const auto consider_radius = geometry_utilities.IsValuePositive(ball_radius,
                                                                        geometry_utilities.Tolerance1D());

        std::vector<double> borders_distance;
        if (consider_radius)
        {
          std::vector<unsigned int> border_sources;
          border_sources.reserve(mesh_graph.mesh_vertices.size());
          for (unsigned int v = 0; v < mesh_graph.mesh_vertices.size(); ++v)
          {
            const auto& vertex = mesh_graph.mesh_vertices.at(v);
            switch (vertex.type)
            {
              case vertex_type::Internal:
                continue;
              case vertex_type::Border:
                border_sources.push_back(v);
                break;
              default:
                throw std::runtime_error("Unknown vertex type");
            }
          }
          border_sources.shrink_to_fit();

          const auto borders_distance_config = Gedim::Algorithms::FastMarchingMethod::create_configuration(geometry_utilities,
                                                                                                           mesh_graph,
                                                                                                           border_sources,
                                                                                                           0.0,
                                                                                                           1.0,
                                                                                                           eikonal_minimization_type);
          const auto borders_distance_computation = Gedim::Algorithms::FastMarchingMethod::fast_marching_method_2D(geometry_utilities,
                                                                                                                   mesh_graph,
                                                                                                                   borders_distance_config);
          borders_distance = borders_distance_computation.vertices_time;
        }

        Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig result;

        result.speed_type = speed_types::Constant;
        result.speed_value = constant_speed;
        result.eikonal_min_type = eikonal_minimization_type;

        result.cell0Ds_initial_condition.resize(mesh_graph.mesh_vertices.size());
        for (unsigned int v = 0; v < mesh_graph.mesh_vertices.size(); ++v)
        {
          auto& initial_condition = result.cell0Ds_initial_condition.at(v);
          initial_condition.type = vertex_initialization_types::Unknown;
          initial_condition.initial_time = infty;
          initial_condition.border_distance = consider_radius ? borders_distance.at(v) : 0.0;
        }

        result.cell1Ds_restriction.resize(mesh_graph.mesh_edges.size());
        for (unsigned int e = 0; e < mesh_graph.mesh_edges.size(); ++e)
        {
          auto& restriction = result.cell1Ds_restriction.at(e);
          restriction.type = edge_restriction_types::Free;
          restriction.reachable_zone = { 0.0, 1.0 };
        }

        std::list<unsigned int> final_sources;
        for (const auto source_index : sources_index)
        {
          auto& initial_condition_s = result.cell0Ds_initial_condition.at(source_index);

          initial_condition_s.type = vertex_initialization_types::Source;
          initial_condition_s.initial_time = 0.0;

          final_sources.push_back(source_index);
        }

        if (consider_radius)
        {
          for (unsigned int v = 0; v < mesh_graph.mesh_vertices.size(); ++v)
          {
            const auto& mesh_vertex_v = mesh_graph.mesh_vertices.at(v);
            auto& initial_condition_v = result.cell0Ds_initial_condition.at(v);

            switch (initial_condition_v.type)
            {
              case vertex_initialization_types::Source:
                continue;
              case vertex_initialization_types::Unreachable:
              case vertex_initialization_types::Unknown:
              case vertex_initialization_types::Initialized:
                break;
              default:
                throw std::runtime_error("Unknown vertex initial condition type");
            }

            bool distance_from_source_lower_R = false;
            for (const auto source_index : sources_index)
            {
              const auto& mesh_vertex_s = mesh_graph.mesh_vertices.at(source_index);

              const auto distance_from_source = geometry_utilities.PointDistance(mesh_vertex_s.coordinate,
                                                                                 mesh_vertex_v.coordinate);

              if (geometry_utilities.IsValueLowerOrEqual(distance_from_source,
                                                         ball_radius,
                                                         geometry_utilities.Tolerance1D()))
              {
                distance_from_source_lower_R = true;
                break;
              }
            }

            if (distance_from_source_lower_R)
            {
              initial_condition_v.type = vertex_initialization_types::Source;
              initial_condition_v.initial_time = 0.0;

              final_sources.push_back(v);
              continue;
            }

            if (geometry_utilities.IsValueLowerOrEqual(borders_distance.at(v),
                                                       ball_radius,
                                                       geometry_utilities.Tolerance1D()))
            {
              initial_condition_v.type = vertex_initialization_types::Unreachable;
              initial_condition_v.initial_time = infty;
            }
          }

          for (const auto& mesh_edge : mesh_graph.mesh_edges)
          {
            auto& restriction_e = result.cell1Ds_restriction.at(mesh_edge.second);

            const auto e_o = std::min(mesh_edge.first.first,
                                      mesh_edge.first.second);
            const auto e_e = std::max(mesh_edge.first.first,
                                      mesh_edge.first.second);

            const auto& initial_condition_e_o = result.cell0Ds_initial_condition.at(e_o);
            const auto& initial_condition_e_e = result.cell0Ds_initial_condition.at(e_e);

            if (initial_condition_e_o.type == vertex_initialization_types::Unreachable &&
                initial_condition_e_e.type == vertex_initialization_types::Unreachable)
            {
              restriction_e.type = edge_restriction_types::Unreachable;
              continue;
            }

            if (initial_condition_e_o.type == vertex_initialization_types::Unreachable ||
                initial_condition_e_e.type == vertex_initialization_types::Unreachable)
            {
              const auto border_distance_o = borders_distance.at(e_o);
              const auto border_distance_e = borders_distance.at(e_e);

              const double lambda_R = (ball_radius - border_distance_o) / (border_distance_e - border_distance_o);

              if (border_distance_e > border_distance_o)
                restriction_e.reachable_zone = { lambda_R, 1.0 };
              else
                restriction_e.reachable_zone = { 0.0, lambda_R };

              restriction_e.type = edge_restriction_types::Constrained;
              continue;
            }
          }
        }

        for (const auto source_index : final_sources)
        {
          const auto& mesh_vertex_s = mesh_graph.mesh_vertices.at(source_index);
          auto& initial_condition_s = result.cell0Ds_initial_condition.at(source_index);

          for (const auto& vn : mesh_vertex_s.adj_list)
          {
            const auto& mesh_vertex_vn = mesh_graph.mesh_vertices.at(vn);
            auto& initial_condition_vn = result.cell0Ds_initial_condition.at(vn);

            switch (initial_condition_vn.type)
            {
              case vertex_initialization_types::Unreachable:
              case vertex_initialization_types::Source:
                continue;
              case vertex_initialization_types::Unknown:
              case vertex_initialization_types::Initialized:
                break;
              default:
                throw std::runtime_error("Unknown vertex initial condition type");
            }

            const auto time_vn = initial_condition_s.initial_time +
                                 compute_time_ij(geometry_utilities,
                                                 constant_speed,
                                                 mesh_vertex_s.coordinate,
                                                 mesh_vertex_vn.coordinate);

            if (time_vn < initial_condition_vn.initial_time)
            {
              initial_condition_vn.type = vertex_initialization_types::Initialized;
              initial_condition_vn.initial_time = time_vn;
            }
          }
        }

        return result;
      }
      // ***************************************************************************
      Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodResult fast_marching_method_2D(const GeometryUtilities& geometry_utilities,
                                                                                              const MeshGraph& mesh_graph,
                                                                                              const FastMarchingMethodConfig& config)
      {
        if (mesh_graph.mesh_vertices.size() != config.cell0Ds_initial_condition.size() ||
            mesh_graph.mesh_edges.size() != config.cell1Ds_restriction.size())
          throw std::runtime_error("Error in input parameters");

        using pq_item = std::pair<double, unsigned int>;
        using status_types = Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodResult::StatusTypes;

        Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodResult result;

        result.vertices_status.resize(mesh_graph.mesh_vertices.size(), status_types::Far);
        result.vertices_time.resize(mesh_graph.mesh_vertices.size(), infty);

        std::list<unsigned int> starting_trial_vertices;
        std::priority_queue<pq_item, std::vector<pq_item>, std::greater<pq_item>> min_pq;
        for (unsigned int v = 0; v < config.cell0Ds_initial_condition.size(); ++v)
        {
          const auto& initial_condition = config.cell0Ds_initial_condition.at(v);

          switch (initial_condition.type)
          {
            case vertex_initialization_types::Source:
            {
              result.vertices_status.at(v) = status_types::Verified;
              result.vertices_time.at(v) = initial_condition.initial_time;
            }
              break;
            case vertex_initialization_types::Initialized:
            {
              result.vertices_status.at(v) = status_types::Trial;
              result.vertices_time.at(v) = initial_condition.initial_time;
              min_pq.push({ result.vertices_time.at(v), v });
            }
              break;
            case vertex_initialization_types::Unreachable:
            {
              result.vertices_status.at(v) = status_types::Verified;
              result.vertices_time.at(v) = infty;
            }
              break;
            case vertex_initialization_types::Unknown:
            {
              result.vertices_status.at(v) = status_types::Far;
              result.vertices_time.at(v) = infty;
            }
              break;
            default:
              throw std::runtime_error("Unknown vertex initial condition type");
          }
        }

        while (!min_pq.empty())
        {
          const auto selected_v = min_pq.top();
          min_pq.pop();

          if (result.vertices_status.at(selected_v.second) ==
              status_types::Verified)
            continue;

          const auto& mesh_vertex_s = mesh_graph.mesh_vertices.at(selected_v.second);
          result.vertices_status.at(selected_v.second) = status_types::Verified;

          for (const auto& adj_i : mesh_vertex_s.adj_list)
          {
            if (result.vertices_status.at(adj_i) ==
                status_types::Verified)
              continue;

            const auto& mesh_vertex_i = mesh_graph.mesh_vertices.at(adj_i);
            auto new_time_i = std::min(result.vertices_time.at(adj_i),
                                       result.vertices_time.at(selected_v.second) +
                                       compute_time_ij(geometry_utilities,
                                                       config.speed_value,
                                                       mesh_vertex_s.coordinate,
                                                       mesh_vertex_i.coordinate));

            for (const auto k : mesh_graph.mesh_triangles.at(std::make_pair(selected_v.second, adj_i)))
            {
              if (result.vertices_status.at(k) !=
                  status_types::Verified)
                continue;

              const auto& edge_restriction = config.cell1Ds_restriction.at(mesh_graph.mesh_edges.at(std::make_pair(selected_v.second,
                                                                                                                   k)));

              const auto time_k = eikonal_triangle(geometry_utilities,
                                                   mesh_graph,
                                                   config.speed_value,
                                                   result.vertices_time,
                                                   selected_v.second,
                                                   k,
                                                   adj_i,
                                                   edge_restriction,
                                                   config.eikonal_min_type);

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
    } // namespace FastMarchingMethod
  } // namespace Algorithms
} // namespace Gedim
