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
#include "BrentMinimization.hpp"
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
        auto min_f =
            [&](double lambda)
        {
          return lambda * vertices_time.at(i) +
              (1.0 - lambda) * vertices_time.at(j) +
              compute_time_ij(geometry_utilities,
                              constant_speed,
                              mesh_graph.mesh_vertices.at(k).coordinate,
                              lambda * mesh_graph.mesh_vertices.at(i).coordinate +
                              (1.0 - lambda) * mesh_graph.mesh_vertices.at(j).coordinate);
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
          case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::EikonalMinimizationTypes::Simple:
          {
            return std::min({
                              min_f(edge_reachable_zone.at(0)),
                              min_f(0.5 * (edge_reachable_zone.at(1) + edge_reachable_zone.at(0))),
                              min_f(edge_reachable_zone.at(1))
                            });
          }
          case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::EikonalMinimizationTypes::Linear:
          {
            if (geometry_utilities.IsValueZero(std::abs(edge_reachable_zone.at(1) - edge_reachable_zone.at(0)),
                                               geometry_utilities.Tolerance1D()))
              return min_f(edge_reachable_zone.at(0));

            const auto min_result = Gedim::Algorithms::Minimization::brent_minimize(
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
                                                                                           const unsigned int& ball_center_vertex_index,
                                                                                           const double& ball_radius,
                                                                                           const double& constant_speed,
                                                                                           const Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::EikonalMinimizationTypes& eikonal_minimization_type)
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
                                                                                                         1.0,
                                                                                                         Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::EikonalMinimizationTypes::Linear);
        const auto borders_distance = Gedim::Algorithms::FastMarchingMethod::fast_marching_method_2D(geometry_utilities,
                                                                                                     mesh_graph,
                                                                                                     1.0,
                                                                                                     borders_distance_config);

        Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig result;
        result.eikonal_min_type = eikonal_minimization_type;

        result.cell0Ds_initial_condition.resize(mesh_graph.mesh_vertices.size());
        for (unsigned int v = 0; v < mesh_graph.mesh_vertices.size(); ++v)
        {
          auto& initial_condition = result.cell0Ds_initial_condition.at(v);
          initial_condition.type = vertex_initialization_types::Unknown;
          initial_condition.initial_time = infty;
        }

        result.cell1Ds_restriction.resize(mesh_graph.mesh_edges.size());
        for (unsigned int e = 0; e < mesh_graph.mesh_edges.size(); ++e)
        {
          auto& restriction = result.cell1Ds_restriction.at(e);
          restriction.type = edge_restriction_types::Free;
          restriction.reachable_zone = { 0.0, 1.0 };
        }

        for (unsigned int v = 0; v < mesh_graph.mesh_vertices.size(); ++v)
        {
          const auto& vertex = mesh_graph.mesh_vertices.at(v);
        }

        return result;
      }
      // ***************************************************************************
      Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig create_configuration(const GeometryUtilities& geometry_utilities,
                                                                                           const MeshGraph& mesh_graph,
                                                                                           const std::vector<unsigned int>& sources_index,
                                                                                           const double& constant_speed,
                                                                                           const Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::EikonalMinimizationTypes& eikonal_minimization_type)
      {
        Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig result;

        result.eikonal_min_type = eikonal_minimization_type;

        result.cell0Ds_initial_condition.resize(mesh_graph.mesh_vertices.size());
        for (unsigned int v = 0; v < mesh_graph.mesh_vertices.size(); ++v)
        {
          auto& initial_condition = result.cell0Ds_initial_condition.at(v);
          initial_condition.type = vertex_initialization_types::Unknown;
          initial_condition.initial_time = infty;
        }

        result.cell1Ds_restriction.resize(mesh_graph.mesh_edges.size());
        for (unsigned int e = 0; e < mesh_graph.mesh_edges.size(); ++e)
        {
          auto& restriction = result.cell1Ds_restriction.at(e);
          restriction.type = edge_restriction_types::Free;
          restriction.reachable_zone = { 0.0, 1.0 };
        }

        for (const auto source_index : sources_index)
        {
          const auto& mesh_vertex_s = mesh_graph.mesh_vertices.at(source_index);
          auto& restriction_s = result.cell0Ds_initial_condition.at(source_index);

          restriction_s.type = vertex_initialization_types::Source;
          restriction_s.initial_time = 0.0;

          for (const auto& vn : mesh_vertex_s.adj_list)
          {
            const auto& mesh_vertex_vn = mesh_graph.mesh_vertices.at(vn);
            auto& restriction_vn = result.cell0Ds_initial_condition.at(vn);

            if (restriction_vn.type ==
                vertex_initialization_types::Source)
              continue;

            const auto time_vn = restriction_s.initial_time +
                                 compute_time_ij(geometry_utilities,
                                                 constant_speed,
                                                 mesh_vertex_s.coordinate,
                                                 mesh_vertex_vn.coordinate);

            if (time_vn < restriction_vn.initial_time)
            {
              restriction_vn.type = vertex_initialization_types::Initialized;
              restriction_vn.initial_time = time_vn;
            }
          }
        }

        return result;
      }
      // ***************************************************************************
      Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodResult fast_marching_method_2D(const GeometryUtilities& geometry_utilities,
                                                                                              const MeshGraph& mesh_graph,
                                                                                              const double& constant_speed,
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
                                                       constant_speed,
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
                                                   constant_speed,
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
