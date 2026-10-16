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

#include "FileTextReader.hpp"
#include "MeshUtilities.hpp"

namespace Gedim
{
// ***************************************************************************
void MeshUtilities::AddHangingNodesToMesh(const std::map<unsigned int, std::vector<double>> &hanging_per_edges, Gedim::IMeshDAO &mesh)
{
    std::map<unsigned int, std::vector<unsigned int>> edge_children;
    std::map<unsigned int, std::vector<unsigned int>> edge_new_vertices;

    for (const auto& edge_pair : hanging_per_edges)
    {
        const unsigned int edge_id = edge_pair.first;
        const Eigen::Vector3d origin = mesh.Cell1DOriginCoordinates(edge_id);
        const Eigen::Vector3d end = mesh.Cell1DEndCoordinates(edge_id);

        const auto &curvilinear = edge_pair.second;
        const unsigned int num_hanging_nodes = curvilinear.size();

        unsigned int cell0DIndex = mesh.Cell0DAppend(num_hanging_nodes);
        unsigned int cell1DIndex = mesh.Cell1DAppend(num_hanging_nodes);

        unsigned int origin_id = mesh.Cell1DOrigin(edge_id);
        unsigned int end_id = mesh.Cell1DEnd(edge_id);
        for (unsigned int p = 0; p < num_hanging_nodes; p++)
        {
            const Eigen::Vector3d coordinate = origin + curvilinear[p] * (end - origin);

            mesh.Cell0DSetState(cell0DIndex, true);
            mesh.Cell0DInsertCoordinates(cell0DIndex, coordinate);
            mesh.Cell0DSetMarker(cell0DIndex, mesh.Cell1DMarker(edge_id));
            cell0DIndex++;

            if (p == 0)
            {
                mesh.Cell1DInsertExtremes(edge_id, origin_id, cell0DIndex - 1);
                origin_id = cell0DIndex - 1;

                edge_children.insert({edge_id, {edge_id}});
                edge_new_vertices.insert({edge_id, {cell0DIndex - 1}});
            }
            else
            {
                mesh.Cell1DSetState(cell1DIndex, true);
                mesh.Cell1DInsertExtremes(cell1DIndex, origin_id, cell0DIndex - 1);
                mesh.Cell1DSetMarker(cell1DIndex, mesh.Cell1DMarker(edge_id));
                cell1DIndex++;
                origin_id = cell0DIndex - 1;

                edge_children[edge_id].push_back(cell1DIndex - 1);
                edge_new_vertices[edge_id].push_back(cell0DIndex - 1);
            }
        }

        mesh.Cell1DSetState(cell1DIndex, true);
        mesh.Cell1DInsertExtremes(cell1DIndex, origin_id, end_id);
        mesh.Cell1DSetMarker(cell1DIndex, mesh.Cell1DMarker(edge_id));
        edge_children[edge_id].push_back(cell1DIndex);
        cell1DIndex++;
    }

    for (unsigned int f = 0; f < mesh.Cell2DTotalNumber(); f++)
    {
        std::vector<unsigned int> vertices;
        std::vector<unsigned int> edges;

        bool modified = false;

        for (unsigned int e = 0; e < mesh.Cell2DNumberEdges(f); e++)
        {
            const unsigned int edge_id = mesh.Cell2DEdge(f, e);
            const auto it = edge_children.find(edge_id);

            if (it == edge_children.end())
            {
                edges.push_back(edge_id);
                vertices.push_back(mesh.Cell2DVertex(f, e));
            }
            else
            {
                modified = true;
                const std::vector<unsigned int> new_vertices = edge_new_vertices[edge_id];
                const std::vector<unsigned int> new_edges = edge_children[edge_id];

                const unsigned int last_vertex = mesh.Cell2DVertex(f, e);
                vertices.push_back(last_vertex);
                const unsigned int old_origin = mesh.Cell1DOrigin(edge_id);

                if (old_origin == last_vertex)
                {
                    for (unsigned int new_edge = 0; new_edge < new_vertices.size(); new_edge++)
                    {
                        vertices.push_back(new_vertices[new_edge]);
                        edges.push_back(new_edges[new_edge]);
                    }

                    edges.push_back(new_edges.back());
                }
                else
                {
                    for (unsigned int new_edge = new_edges.size() - 1; new_edge >= 1; new_edge--)
                    {
                        vertices.push_back(new_vertices[new_edge - 1]);
                        edges.push_back(new_edges[new_edge]);
                    }

                    edges.push_back(new_edges[0]);
                }
            }
        }

        if (modified)
        {
            mesh.Cell2DInitializeVertices(f, vertices.size());
            mesh.Cell2DInitializeEdges(f, edges.size());
            mesh.Cell2DInsertVertices(f, vertices);
            mesh.Cell2DInsertEdges(f, edges);
        }
    }

    for (unsigned int c = 0; c < mesh.Cell3DTotalNumber(); c++)
    {
        std::vector<unsigned int> vertices = mesh.Cell3DVertices(c);
        std::vector<unsigned int> edges;

        bool modified = false;

        for (unsigned int e = 0; e < mesh.Cell3DNumberEdges(c); e++)
        {
            const unsigned int edge_id = mesh.Cell3DEdge(c, e);
            const auto it = edge_children.find(edge_id);

            if (it != edge_children.end())
            {
                modified = true;
                const std::vector<unsigned int> new_vertices = edge_new_vertices[edge_id];
                const std::vector<unsigned int> new_edges = it->second;

                for (unsigned int new_edge = 0; new_edge < new_vertices.size(); new_edge++)
                {
                    vertices.push_back(new_vertices[new_edge]);
                    edges.push_back(new_edges[new_edge]);
                }

                edges.push_back(new_edges.back());
            }
            else
                edges.push_back(edge_id);
        }

        if (modified)
        {
            mesh.Cell3DAddEdges(c, edges);
            mesh.Cell3DAddVertices(c, vertices);
        }
    }
}

// ***************************************************************************
} // namespace Gedim
