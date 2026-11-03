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
      struct MeshGraph final
      {
          struct Vertex
          {
              enum struct VertexTypes
              {
                Internal = 0,
                Border = 1
              };

              Gedim::Algorithms::FastMarchingMethod::MeshGraph::Vertex::VertexTypes type;
              Eigen::Vector3d coordinate;
              std::list<unsigned int> adj_list;
          };

          std::vector<Vertex> mesh_vertices;
          std::map<std::pair<unsigned int, unsigned int>, unsigned int> mesh_edges;
          std::map<std::pair<unsigned int, unsigned int>, std::list<unsigned int>> mesh_triangles;
      };

      struct FastMarchingMethodConfig
      {
          enum struct EikonalMinimizationTypes
          {
            Simple = 0, // only three points are considered
            Linear = 1 // minimization on the line is considered
          };

          struct VertexInitialCondition
          {
              enum struct InitializationTypes
              {
                Source = 0,
                Initialized = 1,
                Unreachable = 2,
                Unknown = 3
              };

              Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::VertexInitialCondition::InitializationTypes type;
              double initial_time;
          };

          struct EdgeRestriction
          {
              enum struct RestrictionTypes
              {
                Free = 0,
                Constrained = 1,
                Unreachable = 2
              };

              Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::EdgeRestriction::RestrictionTypes type;
              std::array<double, 2> reachable_zone;
          };

          Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::EikonalMinimizationTypes eikonal_min_type;
          std::vector<VertexInitialCondition> cell0Ds_initial_condition;
          std::vector<EdgeRestriction> cell1Ds_restriction;
      };

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

      Gedim::Algorithms::FastMarchingMethod::MeshGraph create_mesh_graph(const Gedim::IMeshDAO& mesh);

      Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig create_configuration(const GeometryUtilities& geometry_utilities,
                                                                                           const MeshGraph& mesh_graph,
                                                                                           const std::vector<unsigned int>& sources_index,
                                                                                           const double& constant_speed,
                                                                                           const Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::EikonalMinimizationTypes& eikonal_minimization_type);

      Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig create_configuration(const GeometryUtilities& geometry_utilities,
                                                                                           const MeshGraph& mesh_graph,
                                                                                           const std::vector<unsigned int>& sources_index,
                                                                                           const double& constant_speed,
                                                                                           const Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodConfig::EikonalMinimizationTypes& eikonal_minimization_type);

      Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodResult fast_marching_method_2D(const Gedim::GeometryUtilities& geometry_utilities,
                                                                                              const MeshGraph& mesh_graph,
                                                                                              const double& constant_speed,
                                                                                              const FastMarchingMethodConfig& config);
    } // namespace FastMarchingMethod
  } // namespace Algorithms
} // namespace Gedim

#endif // __Fast_Marching_Method_H
