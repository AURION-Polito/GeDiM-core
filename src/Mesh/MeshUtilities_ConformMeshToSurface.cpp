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
  Gedim::MeshUtilities::ConformMeshToSurfaceResult Gedim::MeshUtilities::ConformMeshToSurface(
      const Gedim::MeshUtilities::IntersectMeshWithSurfaceResult &surface_intersections,
      const Gedim::MeshUtilities::MarchingSquaresResult::SurfaceSegmentsData &surface_segments,
      const Gedim::MeshUtilities::MarchingCubesResult::SurfaceTrianglesData &surface_triangles,
      const Gedim::IMeshDAO &original_mesh,
      const Gedim::MeshUtilities::MeshSurfaceIntersectionGraphData& surface_intersection_graph,
      Gedim::IMeshDAO &conformed_mesh)
  {
    Gedim::MeshUtilities::ConformMeshToSurfaceResult result;

    if (surface_intersection_graph.volumes.empty())
      conformed_mesh.InitializeDimension(2);
    else
      conformed_mesh.InitializeDimension(3);

    if (!surface_intersection_graph.vertices.empty())
    {
      const unsigned int numCell0Ds = surface_intersection_graph.vertices.size();

      result.new_cell0Ds_to_original_cell0D.resize(numCell0Ds,
                                                   Gedim::MeshUtilities::None);
      result.new_cell0Ds_to_surface_vertex.resize(numCell0Ds,
                                                  Gedim::MeshUtilities::None);
      result.original_cell0Ds_to_new_cell0D.resize(original_mesh.Cell0DTotalNumber(),
                                                   Gedim::MeshUtilities::None);

      conformed_mesh.Cell0DsInitialize(numCell0Ds);
      for (unsigned int v = 0; v < numCell0Ds; v++)
      {
        const auto& vertex = surface_intersection_graph.vertices.at(v);

        if (vertex.original_mesh_cell0D_index != Gedim::MeshUtilities::None)
        {
          result.new_cell0Ds_to_original_cell0D.at(v) = vertex.original_mesh_cell0D_index;
          result.original_cell0Ds_to_new_cell0D.at(vertex.original_mesh_cell0D_index) = v;

          conformed_mesh.Cell0DSetMarker(v,
                                         original_mesh.Cell0DMarker(vertex.original_mesh_cell0D_index));
        }
        else
          conformed_mesh.Cell0DSetMarker(v, 0);

        if (vertex.original_surface_point_index != Gedim::MeshUtilities::None)
          result.new_cell0Ds_to_surface_vertex.at(v) = vertex.original_surface_point_index;

        conformed_mesh.Cell0DInsertCoordinates(v, vertex.coordinate);
        conformed_mesh.Cell0DSetState(v, true);
      }
    }

    if (!surface_intersection_graph.edges.empty())
    {
      const unsigned int numCell1Ds = surface_intersection_graph.edges.size();

      result.new_cell1Ds_to_original_cell1D.resize(numCell1Ds,
                                                   Gedim::MeshUtilities::None);
      result.new_cell1Ds_to_surface_segment.resize(numCell1Ds,
                                                   Gedim::MeshUtilities::None);
      result.original_cell1Ds_to_new_cell1Ds.resize(original_mesh.Cell1DTotalNumber());

      conformed_mesh.Cell1DsInitialize(numCell1Ds);
      for (unsigned int e = 0; e < surface_intersection_graph.edges.size(); e++)
      {
        const auto& edge = surface_intersection_graph.edges.at(e);

        if (edge.original_mesh_cell1D_index != Gedim::MeshUtilities::None)
        {
          conformed_mesh.Cell1DSetMarker(e,
                                         original_mesh.Cell1DMarker(edge.original_mesh_cell1D_index));
        }
        else
          conformed_mesh.Cell1DSetMarker(e, 0);

        conformed_mesh.Cell1DInsertExtremes(e, edge.extremes.at(0), edge.extremes.at(1));
        conformed_mesh.Cell1DSetState(e, true);
      }
    }

    if (!surface_intersection_graph.faces.empty())
    {
      const unsigned int numCell2Ds = surface_intersection_graph.faces.size();

      result.new_cell2Ds_to_original_cell2D.resize(numCell2Ds,
                                                   Gedim::MeshUtilities::None);
      result.new_cell2Ds_to_surface_triangle.resize(numCell2Ds,
                                                    Gedim::MeshUtilities::None);
      result.original_cell2Ds_to_new_cell2Ds.resize(original_mesh.Cell2DTotalNumber());

      conformed_mesh.Cell2DsInitialize(numCell2Ds);
      std::vector<unsigned int> cell2Ds_num_vertices(numCell2Ds);

      for (unsigned int f = 0; f < numCell2Ds; f++)
      {
        const auto& face = surface_intersection_graph.faces.at(f);
        cell2Ds_num_vertices.at(f) = face.intersection_graph_vertices.size();
      }

      conformed_mesh.Cell2DsInitializeVertices(cell2Ds_num_vertices);
      conformed_mesh.Cell2DsInitializeEdges(cell2Ds_num_vertices);

      for (unsigned int f = 0; f < numCell2Ds; f++)
      {
        const auto& face = surface_intersection_graph.faces.at(f);

        if (face.original_mesh_cell2D_index != Gedim::MeshUtilities::None)
        {
          conformed_mesh.Cell2DSetMarker(f,
                                         original_mesh.Cell2DMarker(face.original_mesh_cell2D_index));
        }
        else
          conformed_mesh.Cell2DSetMarker(f, 0);

        conformed_mesh.Cell2DInsertVertices(f, face.intersection_graph_vertices);
        conformed_mesh.Cell2DInsertEdges(f, face.intersection_graph_edges);
        conformed_mesh.Cell2DSetState(f, true);
      }
    }

    if (!surface_intersection_graph.volumes.empty())
    {
      const unsigned int numCell3Ds = surface_intersection_graph.volumes.size();

      result.new_cell3Ds_to_original_cell3D.resize(numCell3Ds,
                                                   Gedim::MeshUtilities::None);
      result.original_cell3Ds_to_new_cell3Ds.resize(original_mesh.Cell3DTotalNumber());

      conformed_mesh.Cell3DsInitialize(numCell3Ds);
      std::vector<unsigned int> cell3DsNumVertices(numCell3Ds);
      std::vector<unsigned int> cell3DsNumEdges(numCell3Ds);
      std::vector<unsigned int> cell3DsNumFaces(numCell3Ds);

      for (unsigned int c = 0; c < numCell3Ds; c++)
      {
        const auto volume = surface_intersection_graph.volumes.at(c);

        cell3DsNumVertices.at(c) = volume.intersection_graph_vertices.size();
        cell3DsNumEdges.at(c) = volume.intersection_graph_edges.size();
        cell3DsNumFaces.at(c) = volume.intersection_graph_faces.size();
      }

      conformed_mesh.Cell3DsInitializeVertices(cell3DsNumVertices);
      conformed_mesh.Cell3DsInitializeEdges(cell3DsNumEdges);
      conformed_mesh.Cell3DsInitializeFaces(cell3DsNumFaces);

      for (unsigned int c = 0; c < numCell3Ds; c++)
      {
        const auto volume = surface_intersection_graph.volumes.at(c);

        if (volume.original_mesh_cell3D_index != Gedim::MeshUtilities::None)
        {
          conformed_mesh.Cell3DSetMarker(c,
                                         original_mesh.Cell3DMarker(volume.original_mesh_cell3D_index));
        }
        else
          conformed_mesh.Cell3DSetMarker(c, 0);


        conformed_mesh.Cell3DInsertVertices(c, volume.intersection_graph_vertices);
        conformed_mesh.Cell3DInsertEdges(c, volume.intersection_graph_edges);
        conformed_mesh.Cell3DInsertFaces(c, volume.intersection_graph_faces);
        conformed_mesh.Cell3DSetState(c, true);
      }
    }

    return result;
  }
  // ***************************************************************************
} // namespace Gedim
