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

#include "MEDIT_Utilities.hpp"
#include "MeshDAOExporterToCsv.hpp"
#include "ObjectFileFormatInterface.hpp"
#include "OpenVolumeMeshInterface.hpp"
#include "UCDUtilities.hpp"
#include "WavefrontOBJInterface.hpp"

using namespace std;
using namespace Eigen;

namespace Gedim
{
// ***************************************************************************
MeshUtilities::MeshGeometricData2DConfig::MeshGeometricData2DConfig(bool Cell2DsBoundingBox,
                                                                    bool Cell2DsTriangulations,
                                                                    bool Cell2DsAreas,
                                                                    bool Cell2DsCentroids,
                                                                    bool Cell2DsDiameters,
                                                                    bool Cell2DsEdgeDirections,
                                                                    bool Cell2DsEdgesCentroid,
                                                                    bool Cell2DsEdgeLengths,
                                                                    bool Cell2DsEdgeTangents,
                                                                    bool Cell2DsEdgeNormals,
                                                                    bool Cell2DsChebyshevCenter,
                                                                    bool Cell2DsTriangulationsByChebyshevCenter,
                                                                    bool Cell2DsInRadius)
    : Cell2DsBoundingBox(Cell2DsBoundingBox), Cell2DsTriangulations(Cell2DsTriangulations), Cell2DsAreas(Cell2DsAreas),
      Cell2DsCentroids(Cell2DsCentroids), Cell2DsDiameters(Cell2DsDiameters),
      Cell2DsEdgeDirections(Cell2DsEdgeDirections), Cell2DsEdgesCentroid(Cell2DsEdgesCentroid),
      Cell2DsEdgeLengths(Cell2DsEdgeLengths), Cell2DsEdgeTangents(Cell2DsEdgeTangents),
      Cell2DsEdgeNormals(Cell2DsEdgeNormals), Cell2DsChebyshevCenter(Cell2DsChebyshevCenter),
      Cell2DsTriangulationsByChebyshevCenter(Cell2DsTriangulationsByChebyshevCenter), Cell2DsInRadius(Cell2DsInRadius)
{
}
// ***************************************************************************
void MeshUtilities::ExtractActiveMesh(Gedim::IMeshDAO &mesh, ExtractActiveMeshData &extractionData)
{
    // remove inactive Cell0Ds
    unsigned int numNewCell0Ds = 0;
    list<unsigned int> cell0DIdToRemove;
    for (unsigned int c = 0; c < mesh.Cell0DTotalNumber(); c++)
    {
        if (!mesh.Cell0DIsActive(c))
        {
            cell0DIdToRemove.push_back(c);
            continue;
        }

        extractionData.NewCell0DToOldCell0D.insert(pair<unsigned int, unsigned int>(numNewCell0Ds, c));
        extractionData.OldCell0DToNewCell0D.insert(pair<unsigned int, unsigned int>(c, numNewCell0Ds));
        numNewCell0Ds++;
    }

    unsigned int removedCell0Ds = 0;
    for (const unsigned int &c : cell0DIdToRemove)
    {
        mesh.Cell0DRemove(c - removedCell0Ds);
        removedCell0Ds++;
    }

    // remove inactive Cell1D
    unsigned int numNewCell1Ds = 0;
    list<unsigned int> cell1DIdToRemove;
    for (unsigned int c = 0; c < mesh.Cell1DTotalNumber(); c++)
    {
        if (!mesh.Cell1DIsActive(c))
        {
            cell1DIdToRemove.push_back(c);
            continue;
        }

        extractionData.NewCell1DToOldCell1D.insert(pair<unsigned int, unsigned int>(numNewCell1Ds, c));
        extractionData.OldCell1DToNewCell1D.insert(pair<unsigned int, unsigned int>(c, numNewCell1Ds));
        numNewCell1Ds++;
    }

    unsigned int removedCell1Ds = 0;
    for (const unsigned int &c : cell1DIdToRemove)
    {
        mesh.Cell1DRemove(c - removedCell1Ds);
        removedCell1Ds++;
    }

    // remove inactive Cell2Ds
    unsigned int numNewCell2Ds = 0;
    list<unsigned int> cell2DIdToRemove;
    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
    {
        if (!mesh.Cell2DIsActive(c))
        {
            cell2DIdToRemove.push_back(c);
            continue;
        }

        extractionData.NewCell2DToOldCell2D.insert(pair<unsigned int, unsigned int>(numNewCell2Ds, c));
        extractionData.OldCell2DToNewCell2D.insert(pair<unsigned int, unsigned int>(c, numNewCell2Ds));
        numNewCell2Ds++;
    }

    unsigned int removedCell2Ds = 0;
    for (const unsigned int &c : cell2DIdToRemove)
    {
        mesh.Cell2DRemove(c - removedCell2Ds);
        removedCell2Ds++;
    }

    // remove inactive Cell3Ds
    unsigned int numNewCell3Ds = 0;
    list<unsigned int> cell3DIdToRemove;
    for (unsigned int c = 0; c < mesh.Cell3DTotalNumber(); c++)
    {
        if (!mesh.Cell3DIsActive(c))
        {
            cell3DIdToRemove.push_back(c);
            continue;
        }

        extractionData.NewCell3DToOldCell3D.insert(pair<unsigned int, unsigned int>(numNewCell3Ds, c));
        extractionData.OldCell3DToNewCell3D.insert(pair<unsigned int, unsigned int>(c, numNewCell3Ds));
        numNewCell3Ds++;
    }

    unsigned int removedCell3Ds = 0;
    for (const unsigned int &c : cell3DIdToRemove)
    {
        mesh.Cell3DRemove(c - removedCell3Ds);
        removedCell3Ds++;
    }

    mesh.Compress();
}
// ***************************************************************************
MeshUtilities::FilterMeshData MeshUtilities::FilterActiveMesh(const Gedim::IMeshDAO &mesh)
{
    std::list<unsigned int> active_cell0Ds;
    std::list<unsigned int> active_cell1Ds;
    std::list<unsigned int> active_cell2Ds;
    std::list<unsigned int> active_cell3Ds;

    for (unsigned int c = 0; c < mesh.Cell0DTotalNumber(); c++)
    {
        if (mesh.Cell0DIsActive(c))
            active_cell0Ds.push_back(c);
    }

    for (unsigned int c = 0; c < mesh.Cell1DTotalNumber(); c++)
    {
        if (mesh.Cell1DIsActive(c))
            active_cell1Ds.push_back(c);
    }

    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
    {
        if (mesh.Cell2DIsActive(c))
            active_cell2Ds.push_back(c);
    }

    for (unsigned int c = 0; c < mesh.Cell3DTotalNumber(); c++)
    {
        if (mesh.Cell3DIsActive(c))
            active_cell3Ds.push_back(c);
    }

    return {std::vector<unsigned int>(active_cell0Ds.begin(), active_cell0Ds.end()),
            std::vector<unsigned int>(active_cell1Ds.begin(), active_cell1Ds.end()),
            std::vector<unsigned int>(active_cell2Ds.begin(), active_cell2Ds.end()),
            std::vector<unsigned int>(active_cell3Ds.begin(), active_cell3Ds.end())};
}
// ***************************************************************************
void MeshUtilities::ImportOpenVolumeMesh(const std::string &ovmFilePath, Gedim::IMeshDAO &mesh, std::vector<std::vector<bool>> &meshCell3DsFacesOrientation)
{
    OpenVolumeMeshInterface openVolumeMeshInterface;
    openVolumeMeshInterface.ImportMeshFromFile(ovmFilePath, mesh, meshCell3DsFacesOrientation);
}
// ***************************************************************************
void MeshUtilities::ExportMeshToOpenVolume(const Gedim::IMeshDAO &mesh,
                                           const std::vector<std::vector<bool>> &meshCell3DsFacesOrientation,
                                           const std::string &ovmFilePath)
{
    OpenVolumeMeshInterface openVolumeMeshInterface;
    openVolumeMeshInterface.ExportMeshToFile(mesh, meshCell3DsFacesOrientation, ovmFilePath);
}
// ***************************************************************************
void MeshUtilities::ImportObjectFileFormat(const std::string &offFilePath, Gedim::IMeshDAO &mesh)
{
    ObjectFileFormatInterface objectFileFormatInterface;
    objectFileFormatInterface.ImportMeshFromFile(offFilePath, MeshUtilities(), mesh);
}
// ***************************************************************************
void MeshUtilities::ExportMeshToObjectFileFormat(const Gedim::IMeshDAO &mesh, const std::string &offFilePath)
{
    ObjectFileFormatInterface objectFileFormatInterface;
    objectFileFormatInterface.ExportMeshToFile(mesh, offFilePath);
}
// ***************************************************************************
void MeshUtilities::ImportWavefrontOBJ(const std::string &offFilePath, Gedim::IMeshDAO &mesh)
{
    WavefrontOBJInterface wavefrontOBJInterface;
    wavefrontOBJInterface.ImportMeshFromFile(offFilePath, MeshUtilities(), mesh);
}
// ***************************************************************************
void MeshUtilities::ExportMeshToWavefrontOBJ(const Gedim::IMeshDAO &mesh, const std::string &offFilePath)
{
    WavefrontOBJInterface wavefrontOBJInterface;
    wavefrontOBJInterface.ExportMeshToFile(mesh, offFilePath);
}
// ***************************************************************************
void MeshUtilities::ExportMeshToUCD(const Gedim::IMeshDAO &mesh, const std::string &exportFolder, const std::string &fileName, const bool &separateFile)
{
    string cell0DsFolder = exportFolder;
    string cell1DsFolder = exportFolder;
    string cell2DsFolder = exportFolder;
    string cell3DsFolder = exportFolder;

    if (separateFile)
    {
        Gedim::Output::CreateFolder(exportFolder + "/Cell0Ds");
        Gedim::Output::CreateFolder(exportFolder + "/Cell1Ds");
        Gedim::Output::CreateFolder(exportFolder + "/Cell2Ds");
        Gedim::Output::CreateFolder(exportFolder + "/Cell3Ds");

        cell0DsFolder = exportFolder + "/Cell0Ds";
        cell1DsFolder = exportFolder + "/Cell1Ds";
        cell2DsFolder = exportFolder + "/Cell2Ds";
        cell3DsFolder = exportFolder + "/Cell3Ds";
    }

    // Export Cell0Ds
    if (mesh.Cell0DTotalNumber() > 0)
    {
        vector<double> id(mesh.Cell0DTotalNumber());
        vector<double> marker(mesh.Cell0DTotalNumber());
        vector<double> active(mesh.Cell0DTotalNumber());

        for (unsigned int g = 0; g < mesh.Cell0DTotalNumber(); g++)
        {
            id[g] = g;
            marker[g] = mesh.Cell0DMarker(g);
            active[g] = mesh.Cell0DIsActive(g);
        }

        vector<UCDProperty<double>> properties(3 + mesh.Cell0DNumberDoubleProperties());
        vector<vector<double>> propertyValues(mesh.Cell0DNumberDoubleProperties());

        properties[0] = {"Id", "", static_cast<unsigned int>(id.size()), 1, id.data()};
        properties[1] = {"Marker", "", static_cast<unsigned int>(marker.size()), 1, marker.data()};
        properties[2] = {"Active", "", static_cast<unsigned int>(active.size()), 1, active.data()};

        for (unsigned int p = 0; p < mesh.Cell0DNumberDoubleProperties(); p++)
        {
            propertyValues[p].resize(mesh.Cell0DTotalNumber());
            for (unsigned int g = 0; g < mesh.Cell0DTotalNumber(); g++)
            {
                propertyValues[p][g] = mesh.Cell0DDoublePropertySize(g, p) == 1 ? mesh.Cell0DDoublePropertyValue(g, p, 0) : 0.0;
            }

            properties[3 + p] = {mesh.Cell0DDoublePropertyId(p),
                                 "",
                                 static_cast<unsigned int>(propertyValues[p].size()),
                                 1,
                                 propertyValues[p].data()};
        }

        Gedim::UCDUtilities exporter;
        exporter.ExportPoints(cell0DsFolder + "/" + "Cell0Ds_" + fileName + ".inp", mesh.Cell0DsCoordinates(), properties);
    }

    // Export Cell1Ds
    if (mesh.Cell1DTotalNumber() > 0)
    {
        vector<double> id(mesh.Cell1DTotalNumber());
        vector<double> marker(mesh.Cell1DTotalNumber());
        vector<double> active(mesh.Cell1DTotalNumber());

        for (unsigned int g = 0; g < mesh.Cell1DTotalNumber(); g++)
        {
            id[g] = g;
            marker[g] = mesh.Cell1DMarker(g);
            active[g] = mesh.Cell1DIsActive(g);
        }

        vector<UCDProperty<double>> properties(3 + mesh.Cell1DNumberDoubleProperties());
        vector<vector<double>> propertyValues(mesh.Cell1DNumberDoubleProperties());

        properties[0] = {"Id", "", static_cast<unsigned int>(id.size()), 1, id.data()};
        properties[1] = {"Marker", "", static_cast<unsigned int>(marker.size()), 1, marker.data()};
        properties[2] = {"Active", "", static_cast<unsigned int>(active.size()), 1, active.data()};

        for (unsigned int p = 0; p < mesh.Cell1DNumberDoubleProperties(); p++)
        {
            propertyValues[p].resize(mesh.Cell1DTotalNumber());
            for (unsigned int g = 0; g < mesh.Cell1DTotalNumber(); g++)
            {
                propertyValues[p][g] = mesh.Cell1DDoublePropertySize(g, p) == 1 ? mesh.Cell1DDoublePropertyValue(g, p, 0) : 0.0;
            }

            properties[3 + p] = {mesh.Cell1DDoublePropertyId(p),
                                 "",
                                 static_cast<unsigned int>(propertyValues[p].size()),
                                 1,
                                 propertyValues[p].data()};
        }

        Gedim::UCDUtilities exporter;
        exporter.ExportSegments(cell1DsFolder + "/" + "Cell1Ds_" + fileName + ".inp",
                                mesh.Cell0DsCoordinates(),
                                mesh.Cell1DsExtremes(),
                                {},
                                properties);
    }

    // Export Cell2Ds
    if (mesh.Cell2DTotalNumber() > 0)
    {
        vector<double> id(mesh.Cell2DTotalNumber());
        vector<double> marker(mesh.Cell2DTotalNumber());
        vector<double> active(mesh.Cell2DTotalNumber());

        for (unsigned int g = 0; g < mesh.Cell2DTotalNumber(); g++)
        {
            id[g] = g;
            marker[g] = mesh.Cell2DMarker(g);
            active[g] = mesh.Cell2DIsActive(g);
        }

        vector<UCDProperty<double>> properties(3 + mesh.Cell2DNumberDoubleProperties());
        vector<vector<double>> propertyValues(mesh.Cell2DNumberDoubleProperties());

        properties[0] = {"Id", "", static_cast<unsigned int>(id.size()), 1, id.data()};
        properties[1] = {"Marker", "", static_cast<unsigned int>(marker.size()), 1, marker.data()};
        properties[2] = {"Active", "", static_cast<unsigned int>(active.size()), 1, active.data()};

        for (unsigned int p = 0; p < mesh.Cell2DNumberDoubleProperties(); p++)
        {
            propertyValues[p].resize(mesh.Cell2DTotalNumber());
            for (unsigned int g = 0; g < mesh.Cell2DTotalNumber(); g++)
            {
                propertyValues[p][g] = mesh.Cell2DDoublePropertySize(g, p) == 1 ? mesh.Cell2DDoublePropertyValue(g, p, 0) : 0.0;
            }

            properties[3 + p] = {mesh.Cell2DDoublePropertyId(p),
                                 "",
                                 static_cast<unsigned int>(propertyValues[p].size()),
                                 1,
                                 propertyValues[p].data()};
        }

        Gedim::UCDUtilities exporter;
        exporter.ExportPolygons(cell2DsFolder + "/" + "Cell2Ds_" + fileName + ".inp",
                                mesh.Cell0DsCoordinates(),
                                mesh.Cell2DsVertices(),
                                {},
                                properties);
    }

    // Export Cell3Ds
    if (mesh.Cell3DTotalNumber() > 0)
    {
        vector<double> id(mesh.Cell3DTotalNumber());
        vector<double> marker(mesh.Cell3DTotalNumber());
        vector<double> active(mesh.Cell3DTotalNumber());

        for (unsigned int g = 0; g < mesh.Cell3DTotalNumber(); g++)
        {
            id[g] = g;
            marker[g] = mesh.Cell3DMarker(g);
            active[g] = mesh.Cell3DIsActive(g);
        }

        vector<UCDProperty<double>> properties(3 + mesh.Cell3DNumberDoubleProperties());
        vector<vector<double>> propertyValues(mesh.Cell3DNumberDoubleProperties());

        properties[0] = {"Id", "", static_cast<unsigned int>(id.size()), 1, id.data()};
        properties[1] = {"Marker", "", static_cast<unsigned int>(marker.size()), 1, marker.data()};
        properties[2] = {"Active", "", static_cast<unsigned int>(active.size()), 1, active.data()};

        for (unsigned int p = 0; p < mesh.Cell3DNumberDoubleProperties(); p++)
        {
            propertyValues[p].resize(mesh.Cell3DTotalNumber());
            for (unsigned int g = 0; g < mesh.Cell3DTotalNumber(); g++)
            {
                propertyValues[p][g] = mesh.Cell3DDoublePropertySize(g, p) == 1 ? mesh.Cell3DDoublePropertyValue(g, p, 0) : 0.0;
            }

            properties[3 + p] = {mesh.Cell3DDoublePropertyId(p),
                                 "",
                                 static_cast<unsigned int>(propertyValues[p].size()),
                                 1,
                                 propertyValues[p].data()};
        }

        Gedim::UCDUtilities exporter;
        exporter.ExportPolyhedra(cell3DsFolder + "/" + "Cell3Ds_" + fileName + ".inp",
                                 mesh.Cell0DsCoordinates(),
                                 mesh.Cell3DsVertices(),
                                 {},
                                 properties);
    }
}
// ***************************************************************************
void MeshUtilities::ExportMeshToMEDIT(const Gedim::IMeshDAO &mesh,
                                      const std::string &exportFolder,
                                      const std::string &fileName,
                                      const bool &separateFile,
                                      const unsigned int &order)
{
    if (mesh.Cell0DTotalNumber() == 0)
    {
        return;
    }

    std::string cell0DsFolder = exportFolder;
    std::string cell1DsFolder = exportFolder;
    std::string cell2DsFolder = exportFolder;
    std::string cell3DsFolder = exportFolder;

    if (separateFile)
    {
        Gedim::Output::CreateFolder(exportFolder + "/Cell0Ds");
        Gedim::Output::CreateFolder(exportFolder + "/Cell1Ds");
        Gedim::Output::CreateFolder(exportFolder + "/Cell2Ds");
        Gedim::Output::CreateFolder(exportFolder + "/Cell3Ds");

        cell0DsFolder = exportFolder + "/Cell0Ds";
        cell1DsFolder = exportFolder + "/Cell1Ds";
        cell2DsFolder = exportFolder + "/Cell2Ds";
        cell3DsFolder = exportFolder + "/Cell3Ds";
    }

    const auto points_references_id = mesh.Cell0DsMarker();
    const auto points = mesh.Cell0DsCoordinates();

    // Export Cell0Ds
    {

        Gedim::MEDIT_Utilities exporter;
        exporter.ExportPoints(cell0DsFolder + "/" + "Cell0Ds_" + fileName + ".mesh", points, points_references_id);
    }

    // Export Cell1Ds
    if (mesh.Cell1DTotalNumber() > 0)
    {
        const auto cells_references_id = mesh.Cell1DsMarker();

        Gedim::MEDIT_Utilities exporter;
        exporter.ExportSegments(cell1DsFolder + "/" + "Cell1Ds_" + fileName + ".mesh", points, mesh.Cell1DsExtremes(), points_references_id, cells_references_id);
    }

    // Export Cell2Ds
    if (mesh.Cell2DTotalNumber() > 0)
    {
        const auto cells_references_id = mesh.Cell2DsMarker();

        Gedim::MEDIT_Utilities exporter;
        exporter.ExportPolygons(cell2DsFolder + "/" + "Cell2Ds_" + fileName + ".mesh",
                                points,
                                mesh.Cell2DsVertices(),
                                points_references_id,
                                cells_references_id,
                                order);
    }

    // Export Cell3Ds
    if (mesh.Cell3DTotalNumber() > 0)
    {
        const auto faces_references_id = mesh.Cell2DsMarker();
        const auto cells_references_id = mesh.Cell3DsMarker();

        std::vector<std::vector<unsigned int>> faces_vertices = mesh.Cell2DsVertices();
        std::vector<std::vector<unsigned int>> polyhedra_vertices = mesh.Cell3DsVertices();

        if (order > 1)
        {
            for (unsigned int f = 0; f < mesh.Cell2DTotalNumber(); f++)
            {
                const std::vector<unsigned int> &face_vertices = faces_vertices[f];
                const unsigned int num_vertices = face_vertices.size();
                std::vector<unsigned int> new_physical_face_vertices;
                std::vector<unsigned int> new_hanging_face_vertices;

                bool physical = true;
                bool first_occurence = true;
                int first_physical = 0;
                for (unsigned int v = 0; v < num_vertices; v++)
                {
                    if (mesh.Cell0DNeighbourCell1Ds(face_vertices[v]).size() != 2)
                    {
                        new_physical_face_vertices.push_back(face_vertices[v]);

                        if (first_occurence)
                            first_occurence = false;
                    }
                    else
                    {
                        if (v == 0)
                            physical = false;

                        if (first_occurence)
                            first_physical++;

                        new_hanging_face_vertices.push_back(face_vertices[v]);
                    }
                }

                if (new_hanging_face_vertices.size() % new_physical_face_vertices.size() != 0)
                    throw std::runtime_error("not valid polygon for medit format");

                std::copy(new_physical_face_vertices.begin(), new_physical_face_vertices.end(), faces_vertices[f].begin());
                if (!physical)
                    std::rotate(new_hanging_face_vertices.begin(),
                                new_hanging_face_vertices.begin() + first_physical,
                                new_hanging_face_vertices.end());

                std::copy(new_hanging_face_vertices.begin(),
                          new_hanging_face_vertices.end(),
                          faces_vertices[f].begin() + new_physical_face_vertices.size());
            }

            for (unsigned int c = 0; c < mesh.Cell3DTotalNumber(); c++)
            {
                const std::vector<unsigned int> &cell_vertices = polyhedra_vertices[c];
                const unsigned int num_vertices = cell_vertices.size();
                std::vector<unsigned int> new_physics_cell_vertices;
                std::vector<unsigned int> new_hanging_cell_vertices;

                for (unsigned int v = 0; v < num_vertices; v++)
                {
                    if (mesh.Cell0DNeighbourCell1Ds(cell_vertices[v]).size() != 2)
                        new_physics_cell_vertices.push_back(cell_vertices[v]);
                    else
                        new_hanging_cell_vertices.push_back(cell_vertices[v]);
                }

                std::copy(new_physics_cell_vertices.begin(), new_physics_cell_vertices.end(), polyhedra_vertices[c].begin());
                std::copy(new_hanging_cell_vertices.begin(),
                          new_hanging_cell_vertices.end(),
                          polyhedra_vertices[c].begin() + new_physics_cell_vertices.size());
            }
        }

        Gedim::MEDIT_Utilities exporter;
        exporter.ExportPolyhedra(cell3DsFolder + "/" + "Cell3Ds_" + fileName + ".mesh",
                                 points,
                                 faces_vertices,
                                 polyhedra_vertices,
                                 points_references_id,
                                 faces_references_id,
                                 cells_references_id,
                                 order);
    }
}
// ***************************************************************************
void MeshUtilities::ExportMeshToCsv(const Gedim::IMeshDAO &mesh, const char &separator, const std::string &exportFolderPath)
{
    // export mesh
    Gedim::MeshFromCsvUtilities meshFromCsvUtilities;
    Gedim::MeshFromCsvUtilities::Configuration exportConfiguration;
    exportConfiguration.Separator = separator;
    exportConfiguration.Folder = exportFolderPath;
    Gedim::MeshDAOExporterToCsv exporter(meshFromCsvUtilities);
    exporter.Export(exportConfiguration, mesh);
}
// ***************************************************************************
std::vector<unsigned int> MeshUtilities::MarkCells(const std::function<Eigen::VectorXi(const Eigen::MatrixXd &)> &marking_function,
                                                   const std::vector<Eigen::MatrixXd> &cells_points,
                                                   const unsigned int default_mark)
{
    std::vector<unsigned int> cells_marks(cells_points.size());
    for (unsigned int c = 0; c < cells_points.size(); c++)
    {
        if (cells_points.at(c).size() == 0)
        {
            cells_marks[c] = default_mark;
            continue;
        }

        const auto cell_marks = marking_function(cells_points.at(c));

        cells_marks[c] = cell_marks.maxCoeff();
    }

    return cells_marks;
}
// ***************************************************************************
} // namespace Gedim
