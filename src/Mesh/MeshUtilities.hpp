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

#ifndef __MeshUtilities_H
#define __MeshUtilities_H

#include "GeometryUtilities.hpp"
#include "IMeshDAO.hpp"

namespace Gedim
{
/// \brief MeshUtilities
/// \copyright See top level LICENSE file for details.
struct MeshUtilities
{
    struct CheckMesh2DConfiguration final
    {
        bool Cell0D_CheckCoordinates2D = true;
        bool Cell0D_CheckDuplications = true;
        bool Cell1D_CheckDuplications = true;
        bool Cell1D_CheckNeighbours = true;
        bool Cell1D_CheckMeasure = true;
        bool Cell2D_CheckEdges = true;
        bool Cell2D_CheckDuplications = true;
        bool Cell2D_CheckConvexity = true;
        bool Cell2D_CheckMeasure = true;
    };

    struct CheckMesh3DConfiguration final
    {
        bool Cell0D_CheckDuplications = true;
        bool Cell1D_CheckDuplications = true;
        bool Cell1D_CheckMeasure = true;
        bool Cell2D_CheckEdges = true;
        bool Cell2D_CheckDuplications = true;
        bool Cell2D_CheckConvexity = true;
        bool Cell2D_CheckMeasure = true;
        bool Cell3D_CheckDuplications = true;
        bool Cell3D_CheckEdges = true;
        bool Cell3D_CheckEdgesAreActive = true;
        bool Cell3D_CheckConvexity = true;
        bool Cell3D_CheckMeasure = true;
        bool Cell3D_CheckFaces = true;
    };

    struct CheckMeshGeometricData2DConfiguration final
    {
        bool Cell1D_CheckMeasure = true;
        bool Cell1D_CheckNormals = true;
        bool Cell2D_CheckMeasure = true;
        bool Cell2D_CheckTriangles = true;
        unsigned int Cell1D_QuadratureOrder = 0;
        unsigned int Cell2D_QuadratureOrder = 0;
    };

    struct CheckMeshGeometricData3DConfiguration final
    {
        bool Cell1D_CheckMeasure = true;
        bool Cell1D_CheckNormals = true;
        bool Cell2D_CheckMeasure = true;
        bool Cell2D_CheckTriangles = true;
        bool Cell2D_CheckNormals = true;
        bool Cell3D_CheckMeasure = true;
        bool Cell3D_CheckTetrahedra = true;
        unsigned int Cell1D_QuadratureOrder = 0;
        unsigned int Cell2D_QuadratureOrder = 0;
        unsigned int Cell3D_QuadratureOrder = 0;
    };

    struct ExtractActiveMeshData final
    {
        std::unordered_map<unsigned int, unsigned int> OldCell0DToNewCell0D; ///< each pair is {old Cell0D index, new
        ///< Cell0D index}
        std::unordered_map<unsigned int, unsigned int> OldCell1DToNewCell1D; ///< each pair is {old Cell1D index, new
        ///< Cell1D index}
        std::unordered_map<unsigned int, unsigned int> OldCell2DToNewCell2D; ///< each pair is {old Cell2D index, new
        ///< Cell2D index}
        std::unordered_map<unsigned int, unsigned int> OldCell3DToNewCell3D; ///< each pair is {old Cell3D index, new
        ///< Cell3D index}
        std::unordered_map<unsigned int, unsigned int> NewCell0DToOldCell0D; ///< each pair is {new Cell0D index, old
        ///< Cell0D index}
        std::unordered_map<unsigned int, unsigned int> NewCell1DToOldCell1D; ///< each pair is {new Cell1D index, old
        ///< Cell1D index}
        std::unordered_map<unsigned int, unsigned int> NewCell2DToOldCell2D; ///< each pair is {new Cell2D index, old
        ///< Cell2D index}
        std::unordered_map<unsigned int, unsigned int> NewCell3DToOldCell3D; ///< each pair is {new Cell3D index, old
        ///< Cell3D index}
    };

    struct FilterMeshData final
    {
        std::vector<unsigned int> Cell0Ds;
        std::vector<unsigned int> Cell1Ds;
        std::vector<unsigned int> Cell2Ds;
        std::vector<unsigned int> Cell3Ds;
    };

    struct ExtractMeshData final
    {
        std::vector<unsigned int> OldCell0DToNewCell0D = {}; ///< each element is [old Cell0D index] = new Cell0D index
        std::vector<unsigned int> OldCell1DToNewCell1D = {}; ///< each element is [old Cell1D index] = new Cell1D index
        std::vector<unsigned int> OldCell2DToNewCell2D = {}; ///< each element is [old Cell2D index] = new Cell2D index
        std::vector<unsigned int> OldCell3DToNewCell3D = {}; ///< each element is [old Cell3D index] = new Cell3D index
        std::vector<unsigned int> NewCell0DToOldCell0D = {}; ///< each element is [new Cell0D index] = old Cell0D index
        std::vector<unsigned int> NewCell1DToOldCell1D = {}; ///< each element is [new Cell1D index] = old Cell1D index
        std::vector<unsigned int> NewCell2DToOldCell2D = {}; ///< each element is [new Cell2D index] = old Cell2D index
        std::vector<unsigned int> NewCell3DToOldCell3D = {}; ///< each element is [new Cell3D index] = old Cell3D index
    };

    struct ComputeMesh2DCell1DsResult final
    {
        Eigen::MatrixXi Cell1Ds;              /// Cell1Ds vertices, size 2 x Cell1DTotalNumber()
        std::vector<Eigen::MatrixXi> Cell2Ds; ///< Cell2Ds vertices and edges, size
        ///< Cell2DTotalNumber()x2xCell2DNumberVertices()
    };

    struct ComputeMesh3DAlignedCell1DsResult final
    {

        Eigen::MatrixXi AlignedCell1Ds;
        std::vector<Eigen::MatrixXi> Cell0DsAlignedCell1DsIndex;
        std::vector<Eigen::MatrixXi> Cell1DsAlignedCell1DsIndex;
        std::vector<Eigen::MatrixXi> Cell3DsAlignedCell1DsIndex;
        std::vector<std::vector<unsigned int>> AlignedCell1Ds_SubCell0Ds;
        std::vector<std::vector<unsigned int>> AlignedCell1Ds_SubCell1Ds;
        std::vector<std::vector<unsigned int>> AlignedCell1Ds_Cell3Ds;
    };

    struct MeshGeometricData1D final
    {
        std::vector<Eigen::MatrixXd> Cell1DsBoundingBox;
        std::vector<Eigen::MatrixXd> Cell1DsVertices;  ///< cell1D vertices coordinates
        std::vector<Eigen::Vector3d> Cell1DsTangents;  ///< cell1D tangents
        std::vector<double> Cell1DsLengths;            ///< cell1D lengths
        std::vector<double> Cell1DsSquaredLengths;     ///< cell1D squared lengths
        std::vector<Eigen::Vector3d> Cell1DsCentroids; ///< cell1D centroids
    };

    struct MeshGeometricData2D final
    {
        std::vector<Eigen::MatrixXd> Cell2DsBoundingBox;
        std::vector<Eigen::MatrixXd> Cell2DsVertices;                    ///< cell2D vertices coordinates
        std::vector<std::vector<Eigen::Matrix3d>> Cell2DsTriangulations; ///< cell2D triangulations
        std::vector<double> Cell2DsAreas;                                ///< cell2D areas
        std::vector<Eigen::Vector3d> Cell2DsCentroids;                   ///< cell2D centroids
        std::vector<double> Cell2DsDiameters;                            ///< cell2D diameters
        std::vector<std::vector<bool>> Cell2DsEdgeDirections;            ///< cell2D edge directions
        std::vector<Eigen::MatrixXd> Cell2DsEdgesCentroid;               ///< cell2D edge centroid
        std::vector<Eigen::VectorXd> Cell2DsEdgeLengths;                 ///< cell2D edge lengths
        std::vector<Eigen::MatrixXd> Cell2DsEdgeTangents;                ///< cell2D edge tangents
        std::vector<Eigen::MatrixXd> Cell2DsEdgeNormals;                 ///< cell2D edge normals
        std::vector<Eigen::Vector3d> Cell2DsChebyshevCenter;             ///< cell2D chebycev node
        std::vector<std::vector<Eigen::Matrix3d>> Cell2DsTriangulationsByChebyshevCenter; ///< cell2D triangulations
        std::vector<double> Cell2DsInRadius;                                              ///< cell2D triangulations
    };

    struct MeshGeometricData2DConfig final
    {
        bool Cell2DsBoundingBox = true;
        bool Cell2DsTriangulations = true; ///< cell2D triangulations
        bool Cell2DsAreas = true;          ///< cell2D areas
        bool Cell2DsCentroids = true;      ///< cell2D centroids
        bool Cell2DsDiameters = true;      ///< cell2D diameters
        bool Cell2DsEdgeDirections = true; ///< cell2D edge directions
        bool Cell2DsEdgesCentroid = true;  ///< cell2D edge centroid
        bool Cell2DsEdgeLengths = true;    ///< cell2D edge lengths
        bool Cell2DsEdgeTangents = true;   ///< cell2D edge tangents
        bool Cell2DsEdgeNormals = true;    ///< cell2D edge normals
        bool Cell2DsChebyshevCenter = false;
        bool Cell2DsTriangulationsByChebyshevCenter = false; ///< cell2D triangulations
        bool Cell2DsInRadius = false;                        ///< cell2D triangulations

        MeshGeometricData2DConfig()
        {
            Cell2DsBoundingBox = true;
            Cell2DsTriangulations = true;
            Cell2DsAreas = true;
            Cell2DsCentroids = true;
            Cell2DsDiameters = true;
            Cell2DsEdgeDirections = true;
            Cell2DsEdgesCentroid = true;
            Cell2DsEdgeLengths = true;
            Cell2DsEdgeTangents = true;
            Cell2DsEdgeNormals = true;
            Cell2DsChebyshevCenter = false;
            Cell2DsTriangulationsByChebyshevCenter = false;
            Cell2DsInRadius = false;
        }

        MeshGeometricData2DConfig(bool Cell2DsBoundingBox,
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
                                  bool Cell2DsInRadius);
    };

    struct MeshGeometricData3D final
    {
        std::vector<Eigen::MatrixXd> Cell3DsVertices;
        std::vector<Eigen::MatrixXi> Cell3DsEdges;
        std::vector<std::vector<Eigen::MatrixXi>> Cell3DsFaces;
        std::vector<Eigen::MatrixXd> Cell3DsBoundingBox;
        std::vector<double> Cell3DsVolumes;
        std::vector<double> Cell3DsDiameters;
        std::vector<Eigen::Vector3d> Cell3DsCentroids;
        std::vector<Eigen::MatrixXd> Cell3DsEdgesCentroid;
        std::vector<Eigen::VectorXd> Cell3DsEdgeLengths;
        std::vector<Eigen::MatrixXd> Cell3DsEdgeTangents;
        std::vector<std::vector<bool>> Cell3DsEdgeDirections;
        std::vector<std::vector<Eigen::MatrixXd>> Cell3DsTetrahedronPoints;
        std::vector<std::vector<Eigen::Vector3d>> Cell3DsFacesTranslations;
        std::vector<std::vector<Eigen::Matrix3d>> Cell3DsFacesRotationMatrices;
        std::vector<std::vector<Eigen::Vector3d>> Cell3DsFacesNormals;
        std::vector<std::vector<std::array<Eigen::Vector3d, 2>>> Cell3DsFacesTangents;
        std::vector<std::vector<bool>> Cell3DsFacesNormalDirections;
        std::vector<std::vector<bool>> Cell3DsFacesNormalGlobalDirection;
        std::vector<std::vector<std::array<bool, 2>>> Cell3DsFacesTangentsGlobalDirection;
        std::vector<std::vector<std::vector<bool>>> Cell3DsFacesEdgeDirections;
        std::vector<std::vector<Eigen::MatrixXd>> Cell3DsFaces3DVertices; ///< faces vertices 3D coordinates
        std::vector<std::vector<Eigen::MatrixXd>> Cell3DsFaces2DVertices; ///< faces vertices 2D coordinates
        std::vector<std::vector<std::vector<Eigen::Matrix3d>>> Cell3DsFaces3DTriangulations; ///< faces triangulations
        ///< 2D
        std::vector<std::vector<std::vector<Eigen::Matrix3d>>> Cell3DsFaces2DTriangulations; ///< faces triangulations
        ///< 2D
        std::vector<std::vector<double>> Cell3DsFacesAreas;                   ///< faces areas
        std::vector<std::vector<Eigen::Vector3d>> Cell3DsFaces2DCentroids;    ///< faces centroids
        std::vector<std::vector<double>> Cell3DsFacesDiameters;               ///< faces diameters
        std::vector<std::vector<Eigen::VectorXd>> Cell3DsFacesEdgeLengths;    ///< faces edge lengths
        std::vector<std::vector<Eigen::MatrixXd>> Cell3DsFacesEdge3DTangents; ///< faces edge 3D tangents
        std::vector<std::vector<Eigen::MatrixXd>> Cell3DsFacesEdges3DCentroid;
        std::vector<std::vector<Eigen::MatrixXd>> Cell3DsFacesEdge2DTangents; ///< faces edge 2D tangents
        std::vector<std::vector<Eigen::MatrixXd>> Cell3DsFacesEdges2DCentroid;
        std::vector<std::vector<Eigen::MatrixXd>> Cell3DsFacesEdge2DNormals; ///< faces edge normals
    };

    struct VTPPolyhedron final
    {
        Eigen::MatrixXd Vertices;                               /// size 3xnumVertices
        std::vector<std::vector<unsigned int>> PolyhedronFaces; /// size numFaces x numFaceVertices
    };

    struct AgglomerateTrianglesResult final
    {
        std::vector<unsigned int> RemovedEdges;
        std::vector<unsigned int> VerticesIndex;
        std::vector<unsigned int> EdgesIndex;
    };

    struct AgglomerateMeshFromTriangularMeshResult final
    {
        struct ConcaveCell2D
        {
            unsigned int Cell2DIndex;
            std::vector<unsigned int> ConvexCell2DsIndex;
        };

        std::vector<Gedim::MeshUtilities::AgglomerateMeshFromTriangularMeshResult::ConcaveCell2D> ConcaveCell2Ds;
        std::vector<unsigned int> RemovedCell1Ds;
        std::vector<unsigned int> RemovedCell2Ds;
    };

    struct AgglomerationInformation final
    {
        std::vector<unsigned int> OriginalCell0DToAgglomeratedCell0Ds = {};
        std::vector<unsigned int> OriginalCell1DToAgglomeratedCell1Ds = {};
        std::vector<unsigned int> OriginalCell2DToAgglomeratedCell2Ds = {};
        std::vector<unsigned int> AgglomeratedCell0DToOriginalCell0Ds = {};
        std::vector<std::vector<unsigned int>> AgglomeratedCell1DToOriginalCell1Ds = {};
        std::vector<std::vector<unsigned int>> AgglomeratedCell2DToOriginalCell2Ds = {};
    };

    struct AgglomerateCell1DInformation final
    {
        std::vector<unsigned int> AgglomerateCell1DVertices;
        std::vector<unsigned int> SubCell1DsRemovedVertices;
    };

    struct AgglomerateCell2DInformation final
    {
        std::vector<unsigned int> AgglomerateCell2DVertices;
        std::vector<unsigned int> AgglomerateCell2DEdges;
        std::vector<unsigned int> SubCell2DsRemovedVertices;
        std::vector<unsigned int> SubCell2DsRemovedEdges;
    };

    struct AgglomerateCell3DInformation final
    {
        std::vector<unsigned int> AgglomerateCell3DVertices;
        std::vector<unsigned int> AgglomerateCell3DEdges;
        std::vector<unsigned int> AgglomerateCell3DFaces;
        std::vector<unsigned int> SubCell3DsRemovedVertices;
        std::vector<unsigned int> SubCell3DsRemovedEdges;
        std::vector<unsigned int> SubCell3DsRemovedFaces;
    };

    struct FindConcaveCell3DFacesConvexCell2DResult final
    {
        struct ConvexCell2D final
        {
            unsigned int ConvexCell3DIndex = 0;
            unsigned int ConvexCell3DFaceIndex = 0;
        };

        std::vector<Gedim::MeshUtilities::FindConcaveCell3DFacesConvexCell2DResult::ConvexCell2D> ConcaveCell3DFacesConvexCell2D = {};
    };

    struct Mesh3DPolyhedron final
    {
        std::vector<unsigned int> VerticesIndex;
        std::vector<unsigned int> EdgesIndex;
        std::vector<unsigned int> FacesIndex;
    };

    struct FindPointMeshPositionResult final
    {
        struct PointMeshPosition final
        {
            enum struct Types
            {
                Unknown = 0,
                Outside = 1,
                Cell0D = 2,
                Cell1D = 3,
                Cell2D = 4,
                Cell3D = 5
            };

            Gedim::MeshUtilities::FindPointMeshPositionResult::PointMeshPosition::Types Type;
            unsigned int Cell_index;
        };

        std::vector<Gedim::MeshUtilities::FindPointMeshPositionResult::PointMeshPosition> MeshPositions;
    };

    struct FindPointCell2DResult final
    {
        struct PointCell2DFound final
        {
            unsigned int Cell2D_index;
            Gedim::GeometryUtilities::PointPolygonPositionResult Cell2D_Position;
        };

        std::vector<Gedim::MeshUtilities::FindPointCell2DResult::PointCell2DFound> Cell2Ds_found;
    };

    struct FindPointCell3DResult final
    {
        struct PointCell3DFound final
        {
            unsigned int Cell3D_index;
            Gedim::GeometryUtilities::PointPolyhedronPositionResult Cell3D_Position;
        };

        std::vector<Gedim::MeshUtilities::FindPointCell3DResult::PointCell3DFound> Cell3Ds_found;
    };

    struct Intersect_mesh_polyhedron_result final
    {
        enum struct Types
        {
            None = 1,     ///< No intersection found
            Vertices = 2, ///< Vertices
        };

        struct Polyhedron_Intersection final
        {
            enum struct Types
            {
                Vertex = 1,
                Edge = 2,
                Face = 3,
                Polyhedron = 4
            };

            Gedim::MeshUtilities::Intersect_mesh_polyhedron_result::Polyhedron_Intersection::Types Type;
            unsigned int Geometry_index;
            std::vector<unsigned int> Cell0Ds_index;
            std::vector<unsigned int> Cell1Ds_index;
            std::vector<unsigned int> Cell2Ds_index;
            std::vector<unsigned int> Cell3Ds_index;
        };

        struct Mesh_Intersections final
        {
            std::map<unsigned int, unsigned int> Cell0Ds_intersections;
            std::map<unsigned int, std::vector<unsigned int>> Cell1Ds_intersections;
            std::map<unsigned int, std::vector<unsigned int>> Cell2Ds_intersections;
            std::map<unsigned int, std::vector<unsigned int>> Cell3Ds_intersections;
        };

        Gedim::MeshUtilities::Intersect_mesh_polyhedron_result::Types Type;
        Eigen::MatrixXd Intersections_Coordinates;
        std::vector<Gedim::MeshUtilities::Intersect_mesh_polyhedron_result::Polyhedron_Intersection> Polyhedron_intersections;
        Gedim::MeshUtilities::Intersect_mesh_polyhedron_result::Mesh_Intersections Mesh_intersections;
    };

    struct MarchingCubesData
    {
        std::vector<int> cell0Ds_marked;
        std::vector<int> cell1Ds_new_points;
        Eigen::MatrixXd surface_points;
    };

  public:
    MeshUtilities(){};
    ~MeshUtilities(){};

    /// \brief Extract Active Cells from mesh
    /// \note the resulting mesh has no inactive elements
    static void ExtractActiveMesh(Gedim::IMeshDAO &mesh, Gedim::MeshUtilities::ExtractActiveMeshData &extractionData);

    static Gedim::MeshUtilities::FilterMeshData FilterActiveMesh(const Gedim::IMeshDAO &mesh);

    /// \brief Extract mesh1D cells from a mesh
    static Gedim::MeshUtilities::FilterMeshData FilterMesh1D(const std::vector<unsigned int> &cell1DsFilter,
                                                             const Gedim::IMeshDAO &mesh);

    /// \brief Extract mesh2D cells from a mesh
    static Gedim::MeshUtilities::FilterMeshData FilterMesh2D(const std::vector<unsigned int> &cell2DsFilter,
                                                             const Gedim::IMeshDAO &mesh);

    /// \brief Extract mesh3D cells from a mesh
    static Gedim::MeshUtilities::FilterMeshData FilterMesh3D(const std::vector<unsigned int> &cell3DsFilter,
                                                             const Gedim::IMeshDAO &mesh);

    static Gedim::MeshUtilities::ExtractMeshData ExtractMesh1D(const std::vector<unsigned int> &cell0DsFilter,
                                                               const std::vector<unsigned int> &cell1DsFilter,
                                                               const Gedim::IMeshDAO &originalMesh,
                                                               Gedim::IMeshDAO &mesh);

    static Gedim::MeshUtilities::ExtractMeshData ExtractMesh2D(const std::vector<unsigned int> &cell0DsFilter,
                                                               const std::vector<unsigned int> &cell1DsFilter,
                                                               const std::vector<unsigned int> &cell2DsFilter,
                                                               const Gedim::IMeshDAO &originalMesh,
                                                               Gedim::IMeshDAO &mesh);

    static Gedim::MeshUtilities::ExtractMeshData ExtractMesh3D(const std::vector<unsigned int> &cell0DsFilter,
                                                               const std::vector<unsigned int> &cell1DsFilter,
                                                               const std::vector<unsigned int> &cell2DsFilter,
                                                               const std::vector<unsigned int> &cell3DsFilter,
                                                               const Gedim::IMeshDAO &originalMesh,
                                                               Gedim::IMeshDAO &mesh);

    /// \brief Fill Mesh 1D From segment Coordinates
    /// \param segmentOrigin the segment origin
    /// \param segmentTangent the segment tangent vector
    /// \param coordinates relative coordinates between [0.0, 1.0]
    /// \param mesh the resulting mesh
    static void FillMesh1D(const Gedim::GeometryUtilities &geometryUtilities,
                           const Eigen::Vector3d &segmentOrigin,
                           const Eigen::Vector3d &segmentTangent,
                           const std::vector<double> &coordinates,
                           Gedim::IMeshDAO &mesh);

    /// \brief Fill a Mesh 2D with vertices, edges and polygons
    /// \param cell0Ds the coordinates as Eigen MatrixXd of cell0Ds, size 3xCell0DTotalNumber()
    /// \param cell1Ds the origin and end as Eigen MatrixXd of cell1Ds, size 2xCell1DTotalNumber()
    /// \param cell2Ds the vertices and edges indices of the cell2Ds ordered counterclockwise, size
    /// Cell2DTotalNumber()x2xCell2DNumberVertices()
    static void FillMesh2D(const Eigen::MatrixXd &cell0Ds,
                           const Eigen::MatrixXi &cell1Ds,
                           const std::vector<Eigen::MatrixXi> &cell2Ds,
                           Gedim::IMeshDAO &mesh);

    static void FillMesh3D(const Eigen::MatrixXd &cell0Ds,
                           const Eigen::MatrixXi &cell1Ds,
                           const std::vector<Eigen::MatrixXi> &cell2Ds,
                           const std::vector<Mesh3DPolyhedron> &cell3Ds,
                           Gedim::IMeshDAO &mesh);

    static void FillMesh3D(const Eigen::MatrixXd &cell0Ds,
                           const std::vector<std::pair<unsigned int, unsigned int>> &cell1Ds_vertices,
                           const std::vector<std::vector<unsigned int>> &cell2Ds_vertices,
                           const std::vector<std::vector<unsigned int>> &cell2Ds_edges,
                           const std::vector<std::vector<unsigned int>> &cell3Ds_vertices,
                           const std::vector<std::vector<unsigned int>> &cell3Ds_edges,
                           const std::vector<std::vector<unsigned int>> &cell3Ds_faces,
                           Gedim::IMeshDAO &mesh);

    static void FillMesh3D(const Eigen::MatrixXd &cell0Ds,
                           const std::vector<std::vector<std::vector<unsigned int>>> &cell3Ds_faces_vertices,
                           Gedim::IMeshDAO &mesh);

    static void FixCell2DsOrientation(const Gedim::GeometryUtilities &geometryUtilities, IMeshDAO &mesh);

    static void SetMesh3DMarker(const unsigned int marker, Gedim::IMeshDAO &mesh);

    /// \brief Compute edges in a Mesh 2D with vertices and polygons
    /// \param cell0Ds the coordinates as Eigen MatrixXd of cell0Ds, size 3xCell0DTotalNumber()
    /// \param cell2Ds the vertices indices of the cell2Ds ordered counterclockwise, size
    /// Cell2DTotalNumber()xCell2DNumberVertices() \return the Cell1Ds data
    static Gedim::MeshUtilities::ComputeMesh2DCell1DsResult ComputeMesh2DCell1Ds(const Eigen::MatrixXd &cell0Ds,
                                                                                 const std::vector<Eigen::VectorXi> &cell2Ds);

    /// \brief Check Mesh2D correctness
    /// \param Gedim::GeometryUtilities the geometry utilities
    /// \param convexMesh a convex 2D mesh
    static void CheckMesh2D(const Gedim::MeshUtilities::CheckMesh2DConfiguration &configuration,
                            const Gedim::GeometryUtilities &geometryUtilities,
                            const Gedim::IMeshDAO &convexMesh);

    /// \brief Check Mesh3D correctness
    /// \param Gedim::GeometryUtilities the geometry utilities
    /// \param mesh a 3D mesh
    static void CheckMesh3D(const Gedim::MeshUtilities::CheckMesh3DConfiguration &configuration,
                            const Gedim::GeometryUtilities &geometryUtilities,
                            const Gedim::IMeshDAO &mesh);

    /// \brief Compute edges in a Mesh 2D with vertices and polygons
    /// \param cell0Ds the coordinates as Eigen MatrixXd of cell0Ds, size 3xCell0DTotalNumber()
    /// \param cell2Ds the vertices indices of the cell2Ds ordered counterclockwise, size
    /// Cell2DTotalNumber()xCell2DNumberVertices() \return the Cell1Ds data
    static Gedim::MeshUtilities::ComputeMesh3DAlignedCell1DsResult ComputeMesh3DAlignedCell1Ds(
        const std::vector<std::vector<std::vector<unsigned int>>> &cell3DsAlignedEdgesVertices,
        const std::vector<std::vector<std::vector<unsigned int>>> &cell3DsAlignedEdgesEdges,
        const Gedim::IMeshDAO &mesh);

    /// \brief Check MeshGeometricData3D correctness
    /// \param Gedim::GeometryUtilities the geometry utilities
    /// \param mesh the 3D mesh
    /// \param geometricData the mesh geometric data
    static void CheckMeshGeometricData3D(const Gedim::MeshUtilities::CheckMeshGeometricData3DConfiguration &configuration,
                                         const Gedim::GeometryUtilities &geometryUtilities,
                                         const Gedim::IMeshDAO &mesh,
                                         const Gedim::MeshUtilities::MeshGeometricData3D &geometricData);

    /// \brief Create a Mesh 1D with a segment
    /// \param segmentVertices the segment coordinates, size 3x2
    /// \param vertexMarkers mesh markers of vertices, size 1xNumPolygonVertices()
    static void Mesh1DFromSegment(const Gedim::GeometryUtilities &geometryUtilities,
                                  const Eigen::MatrixXd &segmentVertices,
                                  const std::vector<unsigned int> vertexMarkers,
                                  Gedim::IMeshDAO &mesh);

    /// \brief Create a Mesh 2D with a polygon
    /// \param polygonVertices the polygon coordinates, size 3xNumPolygonVertices()
    /// \param vertexMarkers mesh markers of vertices, size 1xNumPolygonVertices()
    /// \param edgeMarkers mesh markers of edges, size 1xNumPolygonVertices()
    static void Mesh2DFromPolygon(const Eigen::MatrixXd &polygonVertices,
                                  const std::vector<unsigned int> vertexMarkers,
                                  const std::vector<unsigned int> edgeMarkers,
                                  Gedim::IMeshDAO &mesh);

    /// \brief Set the marker on all the mesh 2D elements laying on the line
    /// \param Gedim::GeometryUtilities the geometry utilities
    /// \param lineTangent the line tangent
    /// \param lineOrigin the line origin
    /// \param lineTangentSquaredLength the line tangent squared length
    /// \param marker the marker
    /// \param mesh the mesh
    static void SetMeshMarkersOnLine(const Gedim::GeometryUtilities &geometryUtilities,
                                     const Eigen::Vector3d &lineOrigin,
                                     const Eigen::Vector3d &lineTangent,
                                     const double &lineTangentSquaredLength,
                                     const unsigned int &marker,
                                     Gedim::IMeshDAO &mesh);

    static void SetMeshMarkersOnSegment(const Gedim::GeometryUtilities &geometryUtilities,
                                        const Eigen::Vector3d &segment_origin,
                                        const Eigen::Vector3d &segment_tangent,
                                        const double &segment_tangent_squared_length,
                                        const unsigned int &marker,
                                        Gedim::IMeshDAO &mesh);

    /// \brief Create a Mesh 3D with a polyhedron
    /// \param polyhedronVertices the polyhedron vertices, size 3 x numVertices
    /// \param polyhedronEdges the polyhedron edges, size 2 x numEdges
    /// \param polyhedronFaces the polyhedron face vertices and edges, size numFaces x 2 x numVertices
    /// \param vertexMarkers mesh markers of vertices, size 1xnumVertices
    /// \param edgeMarkers mesh markers of edges, size 1xnumEdges
    /// \param faceMarkers mesh markers of faces, size 1xnumFaces
    static void Mesh3DFromPolyhedron(const Eigen::MatrixXd &polyhedronVertices,
                                     const Eigen::MatrixXi &polyhedronEdges,
                                     const std::vector<Eigen::MatrixXi> &polyhedronFaces,
                                     const std::vector<unsigned int> vertexMarkers,
                                     const std::vector<unsigned int> edgeMarkers,
                                     const std::vector<unsigned int> faceMarkers,
                                     Gedim::IMeshDAO &mesh);

    /// \brief Set the marker on all the mesh 3D elements laying on the plane
    /// \param Gedim::GeometryUtilities the geometry utilities
    /// \param planeNormal the plane normal
    /// \param planeOrigin the plane origin
    /// \param marker the marker
    /// \param mesh the mesh
    static void SetMeshMarkersOnPlane(const Gedim::GeometryUtilities &geometryUtilities,
                                      const Eigen::Vector3d &planeNormal,
                                      const Eigen::Vector3d &planeOrigin,
                                      const unsigned int &marker,
                                      Gedim::IMeshDAO &mesh);

    static void SetMeshMarkersByFaceNormal(const Gedim::GeometryUtilities &geometryUtilities,
                                           const Eigen::Vector3d &normal,
                                           const std::vector<Eigen::Vector3d> &cell2Ds_normal,
                                           const unsigned int &marker,
                                           Gedim::IMeshDAO &mesh);

    /// \warning Only for convex points
    static void SetMeshMarkersOnPolygon(const Gedim::GeometryUtilities &geometryUtilities,
                                        const Eigen::Vector3d &polygon_plane_normal,
                                        const Eigen::Vector3d &polygon_plane_origin,
                                        const Eigen::MatrixXd &polygon_vertices_2D,
                                        const Eigen::Vector3d &polygon_translation,
                                        const Eigen::Matrix3d &polygon_rotation_matrix,
                                        const unsigned int &marker,
                                        Gedim::IMeshDAO &mesh);

    static void SetMeshMarkersOnPolygon(const Gedim::GeometryUtilities &geometryUtilities,
                                        const Eigen::Vector3d &polygon_plane_normal,
                                        const Eigen::Vector3d &polygon_plane_origin,
                                        const Eigen::MatrixXd &polygon_vertices_2D,
                                        const Eigen::Vector3d &polygon_translation,
                                        const Eigen::Matrix3d &polygon_rotation_matrix,
                                        const std::vector<Eigen::Vector3d> &cell1Ds_centroid,
                                        const std::vector<Eigen::Vector3d> &cell2Ds_centroid,
                                        const unsigned int &marker,
                                        Gedim::IMeshDAO &mesh);

    /// \brief Extract the mesh Cell2D Roots
    /// \param mesh the mesh
    /// \return the root cell for each cell2D, size 1xCell2DTotalNumber()
    static std::vector<unsigned int> MeshCell2DRoots(const Gedim::IMeshDAO &mesh);

    /// \brief Fill Mesh1D Geometric Data given a mesh with convex mesh cells
    /// \param convexMesh the convex mesh
    /// \return the MeshGeometricData computed
    static MeshGeometricData1D FillMesh1DGeometricData(const Gedim::GeometryUtilities &geometryUtilities,
                                                       const Gedim::IMeshDAO &convexMesh);

    static MeshGeometricData1D ImportMeshGeometricData1DFromTxt(const std::string &file_path);
    static void ExportMeshGeometricData1DToTxt(const MeshGeometricData1D &mesh_geometric_data, const std::string &file_path);

    /// \brief Fill Mesh2D Geometric Data given a mesh with convex mesh cells
    /// \param convexMesh the convex mesh
    /// \return the MeshGeometricData computed
    static MeshGeometricData2D FillMesh2DGeometricData(const Gedim::GeometryUtilities &geometryUtilities,
                                                       const Gedim::IMeshDAO &convexMesh);

    /// \brief Fill Mesh2D Geometric Data given a mesh with mesh cells type
    /// \param mesh the mesh
    /// \param meshCell2DsPolygonType the cell2D polygon type
    /// \return the MeshGeometricData computed
    static MeshUtilities::MeshGeometricData2D FillMesh2DGeometricData(
        const Gedim::GeometryUtilities &geometryUtilities,
        const Gedim::IMeshDAO &mesh,
        const std::vector<Gedim::GeometryUtilities::PolygonTypes> &meshCell2DsPolygonType,
        const MeshGeometricData2DConfig &config = MeshGeometricData2DConfig());

    /// \brief Fill Mesh2D Geometric Data starting given a mesh with non convex mesh cells and its convex sub-mesh cells
    /// \param mesh the mesh
    /// \param convexMesh the convex mesh cells of mesh
    /// \param meshCell2DToConvexCell2DIndices the collection of convex cell2Ds for each mesh cell2D
    /// \return the MeshGeometricData computed
    static MeshGeometricData2D FillMesh2DGeometricData(const Gedim::GeometryUtilities &geometryUtilities,
                                                       const Gedim::IMeshDAO &mesh,
                                                       const Gedim::IMeshDAO &convexMesh,
                                                       const std::vector<std::vector<unsigned int>> &meshCell2DToConvexCell2DIndices);

    static MeshGeometricData2D ImportMeshGeometricData2DFromTxt(const std::string &file_path);
    static void ExportMeshGeometricData2DToTxt(const MeshGeometricData2D &mesh_geometric_data, const std::string &file_path);

    /// \brief Fill Mesh3D Geometric Data given a mesh with convex mesh cells
    /// \param convexMesh the convex mesh
    /// \return the MeshGeometricData computed
    static MeshGeometricData3D FillMesh3DGeometricData(const Gedim::GeometryUtilities &geometryUtilities,
                                                       const Gedim::IMeshDAO &convexMesh);

    /// \brief Fill Mesh3D Geometric Data starting given a mesh with non convex mesh cells and its convex sub-mesh cells
    /// \param mesh the mesh
    /// \param convexMesh the convex mesh
    /// \param meshCell2DToConvexCell2DIndices the collection of convex cell2Ds for each mesh cell2D
    /// \param meshCell3DToConvexCell3DIndices the collection of convex cell3Ds for each mesh cell3D
    /// \return the MeshGeometricData computed
    static MeshGeometricData3D FillMesh3DGeometricData(const Gedim::GeometryUtilities &geometryUtilities,
                                                       const Gedim::IMeshDAO &mesh,
                                                       const Gedim::IMeshDAO &convexMesh,
                                                       const std::vector<std::vector<unsigned int>> &meshCell3DToConvexCell3DIndices);

    static MeshGeometricData3D FillMesh3DGeometricData(const Gedim::GeometryUtilities &geometryUtilities,
                                                       const Gedim::IMeshDAO &mesh,
                                                       const std::vector<std::vector<Eigen::MatrixXd>> &cell3Ds_tetra_vertices,
                                                       const std::vector<std::vector<Eigen::Matrix3d>> &cell2Ds_triangles_3D_vertices);

    static MeshGeometricData3D ImportMeshGeometricData3DFromTxt(const std::string &file_path);
    static void ExportMeshGeometricData3DToTxt(const MeshGeometricData3D &mesh_geometric_data, const std::string &file_path);

    static void ComputeCell0DCell1DNeighbours(Gedim::IMeshDAO &mesh);
    static void ComputeCell0DCell2DNeighbours(Gedim::IMeshDAO &mesh);
    static void ComputeCell0DCell3DNeighbours(Gedim::IMeshDAO &mesh);

    /// \brief Compute Cell1D Cell2DNeighbours with given mesh data
    /// \param mesh the resulting mesh
    static void ComputeCell1DCell2DNeighbours(Gedim::IMeshDAO &mesh);

    /// \brief Compute Cell1D Cell3DNeighbours with given mesh data
    /// \param mesh the resulting mesh
    static void ComputeCell1DCell3DNeighbours(Gedim::IMeshDAO &mesh);

    /// \brief Compute Cell2D Cell3DNeighbours with given mesh data
    /// \param mesh the resulting mesh
    static void ComputeCell2DCell3DNeighbours(Gedim::IMeshDAO &mesh);

    /// \brief Crete rectange Mesh on rectangle base x height
    /// \param rectangleOrigin the rectangle origin point
    /// \param rectangleBaseTangent the rectangle base tangent vector
    /// \param rectangleHeightTangent the rectangle height tangent vector
    /// \param baseMeshCurvilinearCoordinates the base mesh 1D curvilinear coordinates
    /// \param heightMeshCurvilinearCoordinates the height mesh 1D curvilinear coordinates
    /// \note markers on border are set as { 1, 2, 3, 4, ..., numVertices } for cell0Ds and { 5, 6, 7, 8, ..., 2 *
    /// numVertices } for cell1Ds
    static void CreateRectangleMesh(const Eigen::Vector3d &rectangleOrigin,
                                    const Eigen::Vector3d &rectangleBaseTangent,
                                    const Eigen::Vector3d &rectangleHeightTangent,
                                    const std::vector<double> &baseMeshCurvilinearCoordinates,
                                    const std::vector<double> &heightMeshCurvilinearCoordinates,
                                    Gedim::IMeshDAO &mesh);

    static void CreateStructuredTriangularMesh(const Eigen::Vector3d &rectangleOrigin,
                                               const Eigen::Vector3d &rectangleBaseTangent,
                                               const Eigen::Vector3d &rectangleHeightTangent,
                                               const std::vector<double> &baseMeshCurvilinearCoordinates,
                                               const std::vector<double> &heightMeshCurvilinearCoordinates,
                                               Gedim::IMeshDAO &mesh);

    static void CreateParallelepipedMesh(const Eigen::Vector3d &rectangleOrigin,
                                         const Eigen::Vector3d &rectangleLengthTangent,
                                         const Eigen::Vector3d &rectangleHeightTangent,
                                         const Eigen::Vector3d &rectangleWidthTangent,
                                         const std::vector<double> &lengthMeshCurvilinearCoordinates,
                                         const std::vector<double> &heightMeshCurvilinearCoordinates,
                                         const std::vector<double> &widthMeshCurvilinearCoordinates,
                                         Gedim::IMeshDAO &mesh);

    static void CreateTrianglePlusHangingNodesMesh(const Eigen::Vector3d &rectangleOrigin,
                                                   const Eigen::Vector3d &rectangleBaseTangent,
                                                   const Eigen::Vector3d &rectangleHeightTangent,
                                                   const std::vector<double> &baseMeshCurvilinearCoordinates,
                                                   const std::vector<double> &heightMeshCurvilinearCoordinates,
                                                   const std::vector<unsigned int> &numberOfAddedVerticesForEachRectangle,
                                                   const Gedim::GeometryUtilities &geometryUtilities,
                                                   Gedim::IMeshDAO &mesh);

    static void CreateRectanglePlusHangingNodesMesh(const Eigen::Vector3d &rectangleOrigin,
                                                    const Eigen::Vector3d &rectangleBaseTangent,
                                                    const Eigen::Vector3d &rectangleHeightTangent,
                                                    const std::vector<double> &baseMeshCurvilinearCoordinates,
                                                    const std::vector<double> &heightMeshCurvilinearCoordinates,
                                                    const std::vector<unsigned int> &numberOfAddedVerticesForEachRectangle,
                                                    const Gedim::GeometryUtilities &geometryUtilities,
                                                    Gedim::IMeshDAO &mesh);

    /// \brief Create a quadrilateral mesh starting from a triangular mesh: the mesh is created by dividing triangles
    /// intro 3 quadrilaterals obtained by connecting the centroid to edge midpoints
    static void CreateQuadrilateralMeshFromTriangularMesh(const Gedim::GeometryUtilities &geometryUtilities,
                                                          const Gedim::IMeshDAO &triangular_mesh,
                                                          Gedim::IMeshDAO &mesh);

    static void MakeMeshTriangularFaces(const std::vector<std::vector<unsigned int>> &faces_triangulation, Gedim::IMeshDAO &mesh);

    /// \brief Import triangular 2D mesh (non standard simple format)
    /// \param cell0Ds_file_path cell0Ds coordinates
    /// \param cell2Ds_file_path cell2Ds connectivity
    /// \param marker_file_path markers file
    static void ImportTriangularMesh(const Gedim::GeometryUtilities &geometry_utilities,
                                     const std::string &cell0Ds_file_path,
                                     const std::string &cell2Ds_file_path,
                                     const std::string &marker_file_path,
                                     const char separator,
                                     Gedim::IMeshDAO &mesh);

    /// \brief Import regn_face 3D mesh (non standard polyhedral format)
    /// \param node_file_path cell0Ds coordinates
    /// \param ele_file_path cell2Ds-cell3Ds connectivity
    static void ImportRegnFaceMesh(const Gedim::GeometryUtilities &geometry_utilities,
                                   const std::string &node_file_path,
                                   const std::string &ele_file_path,
                                   IMeshDAO &mesh);

    /// \brief Import 3D mesh from OVM file
    static void ImportOpenVolumeMesh(const std::string &ovmFilePath,
                                     Gedim::IMeshDAO &mesh,
                                     std::vector<std::vector<bool>> &meshCell3DsFacesOrientation);

    /// \brief Export 3D mesh to OVM file
    static void ExportMeshToOpenVolume(const Gedim::IMeshDAO &mesh,
                                       const std::vector<std::vector<bool>> &meshCell3DsFacesOrientation,
                                       const std::string &ovmFilePath);

    /// \brief Import 2D mesh from OFF file
    static void ImportObjectFileFormat(const std::string &offFilePath, Gedim::IMeshDAO &mesh);

    /// \brief Export 2D mesh to OFF file
    static void ExportMeshToObjectFileFormat(const Gedim::IMeshDAO &mesh, const std::string &offFilePath);

    /// \brief Import 2D mesh from OBJ file
    static void ImportWavefrontOBJ(const std::string &offFilePath, Gedim::IMeshDAO &mesh);

    /// \brief Export 2D mesh to OBJ file
    static void ExportMeshToWavefrontOBJ(const Gedim::IMeshDAO &mesh, const std::string &offFilePath);

    /// \brief Change Polygon Mesh Markers from { 1, 2, 3, 4, ..., numVertices } for cell0Ds and { 5, 6, 7, 8, ..., 2 *
    /// numVertices } for cell1Ds to cell0DMarkers and cell1DMarkers \param polygonVertices the 2D polygon vertices,
    /// size 3xnumVertices \param cell0DMarkers the new cell0D markers, size 1xnumPolygonVertices \param cell1DMarkers
    /// the new cell1D markers, size 1xnumPolygonVertices \param mesh the mesh
    static void ChangePolygonMeshMarkers(const Eigen::MatrixXd &polygonVertices,
                                         const std::vector<unsigned int> &cell0DMarkers,
                                         const std::vector<unsigned int> &cell1DMarkers,
                                         Gedim::IMeshDAO &mesh);

    static void ChangePolyhedronMeshMarkers(const Gedim::GeometryUtilities &geometryUtilities,
                                            const Eigen::MatrixXd &polyhedron_vertices,
                                            const Eigen::MatrixXi &polyhedron_edges,
                                            const std::vector<Eigen::MatrixXi> &polyhedron_faces,
                                            const Eigen::MatrixXd &polyhedron_edges_tangent,
                                            const Eigen::VectorXd &polyhedron_edges_length,
                                            const std::vector<Eigen::Vector3d> &polyhedron_faces_normal,
                                            const std::vector<Eigen::MatrixXd> &polyhedron_faces_vertices,
                                            const std::vector<Eigen::MatrixXd> &polyhedron_faces_vertices_2D,
                                            const std::vector<Eigen::Vector3d> &polyhedron_faces_translation,
                                            const std::vector<Eigen::Matrix3d> &polyhedron_faces_rotation_matrix,
                                            const std::vector<unsigned int> &polyhedron_vertices_marker,
                                            const std::vector<unsigned int> &polyhedron_edges_marker,
                                            const std::vector<unsigned int> &polyhedron_faces_marker,
                                            const std::vector<Eigen::Vector3d> &cell1Ds_centroid,
                                            const std::vector<Eigen::Vector3d> &cell2Ds_centroid,
                                            Gedim::IMeshDAO &mesh);

    /// \brief Export Mesh To UCD
    /// \param mesh the mesh
    /// \param exportFolder the folder in which the mesh is exported
    static void ExportMeshToUCD(const Gedim::IMeshDAO &mesh,
                                const std::string &exportFolder,
                                const std::string &fileName,
                                const bool &separateFile = false);

    /// \brief Export Mesh To MEDIT
    /// \param mesh the mesh
    /// \param exportFolder the folder in which the mesh is exported
    /// \param order = 1 or order 2 (mid points added)
    static void ExportMeshToMEDIT(const Gedim::IMeshDAO &mesh,
                                  const std::string &exportFolder,
                                  const std::string &fileName,
                                  const bool &separateFile = false,
                                  const unsigned int &order = 1);

    /// \brief Convert a mesh cell3D to a geometric polydheron
    /// \param mesh a mesh
    /// \param cell3DIndex the cell3D index
    /// \return polyhedron from mesh 3D cell
    static Gedim::GeometryUtilities::Polyhedron MeshCell3DToPolyhedron(const Gedim::IMeshDAO &mesh, const unsigned int &cell3DIndex);

    /// \brief Convert a mesh cell3D to a VTP polydheron
    /// \param mesh a mesh
    /// \param cell3DIndex the cell3D index
    /// \return VTP polyhedron from mesh 3D cell
    static MeshUtilities::VTPPolyhedron MeshCell3DToVTPPolyhedron(const Gedim::IMeshDAO &mesh, const unsigned int &cell3DIndex);

    /// \brief Split cell2D into subcells
    /// \param cell1DIndex the index of Cell1D from 0 to Cell1DTotalNumber()
    /// \param subCell1Ds the list of sub-cells 1D mesh vertices indices, size 2 x numSubCells)
    /// \param mesh the mesh to update
    /// \return the list of new cell1Ds indices, from 0 to Cell1DTotalNumber()
    static std::vector<unsigned int> SplitCell1D(const unsigned int &cell1DIndex, const Eigen::MatrixXi subCell1Ds, Gedim::IMeshDAO &mesh);

    /// \brief Split cell2D into subcells
    /// \param cell2DIndex the index of Cell2D from 0 to Cell2DTotalNumber()
    /// \param subCell2Ds the list of sub-cells 2D mesh vertices and edges indices, size numSubCells x (2 x numVertices)
    /// \param mesh the mesh to update
    /// \return the list of new cell2Ds indices, from 0 to Cell2DTotalNumber()
    static std::vector<unsigned int> SplitCell2D(const unsigned int &cell2DIndex,
                                                 const std::vector<Eigen::MatrixXi> &subCell2Ds,
                                                 Gedim::IMeshDAO &mesh);

    /// \brief Split cell3D into subcells
    /// \param cell3DIndex the index of Cell3D from 0 to Cell3DTotalNumber()
    /// \param subCell3Ds the list of sub-cells 3D mesh vertices and edges indices, size numSubCells x (2 x numVertices)
    /// \param mesh the mesh to update
    /// \return the list of new cell3Ds indices, from 0 to Cell3DTotalNumber()
    static std::vector<unsigned int> SplitCell3D(const unsigned int &cell3DIndex,
                                                 const std::vector<std::vector<unsigned int>> &subCell3DsVertices,
                                                 const std::vector<std::vector<unsigned int>> &subCell3DsEdges,
                                                 const std::vector<std::vector<unsigned int>> &subCell3DsFaces,
                                                 Gedim::IMeshDAO &mesh);

    static Gedim::MeshUtilities::AgglomerateCell1DInformation AgglomerateCell1Ds(const std::unordered_set<unsigned int> &cell1DsIndex,
                                                                                 const Gedim::IMeshDAO &mesh);

    static Gedim::MeshUtilities::AgglomerateCell2DInformation AgglomerateCell2Ds(const Gedim::GeometryUtilities &geometryUtilities,
                                                                                 const std::unordered_set<unsigned int> &cell2DsIndex,
                                                                                 const Gedim::IMeshDAO &mesh);

    static Gedim::MeshUtilities::AgglomerateCell3DInformation AgglomerateCell3Ds(const std::unordered_set<unsigned int> &cell3DsIndex,
                                                                                 const Gedim::IMeshDAO &mesh);

    static unsigned int AgglomerateCell1Ds(const std::unordered_set<unsigned int> &subCell1DsIndex,
                                           const std::vector<unsigned int> &agglomerateCell1DVertices,
                                           const std::vector<unsigned int> &subCell1DsRemovedCell0Ds,
                                           Gedim::IMeshDAO &mesh,
                                           std::vector<std::vector<unsigned int>> &meshCell1DsOriginalCell1Ds,
                                           const bool mantain_neigh2D_order = false);

    static unsigned int AgglomerateCell2Ds(const std::unordered_set<unsigned int> &subCell2DsIndex,
                                           const std::vector<unsigned int> &agglomerateCell2DVertices,
                                           const std::vector<unsigned int> &agglomerateCell2DEdges,
                                           const std::vector<unsigned int> &subCell2DsRemovedCell0Ds,
                                           const std::vector<unsigned int> &subCell2DsRemovedCell1Ds,
                                           Gedim::IMeshDAO &mesh,
                                           std::vector<std::vector<unsigned int>> &meshCell2DsOriginalCell2Ds);

    static unsigned int AgglomerateCell3Ds(const std::unordered_set<unsigned int> &subCell3DsIndex,
                                           const std::vector<unsigned int> &agglomerateCell3DVertices,
                                           const std::vector<unsigned int> &agglomerateCell3DEdges,
                                           const std::vector<unsigned int> &agglomerateCell3DFaces,
                                           const std::vector<unsigned int> &subCell3DsRemovedCell0Ds,
                                           const std::vector<unsigned int> &subCell3DsRemovedCell1Ds,
                                           const std::vector<unsigned int> &subCell3DsRemovedCell2Ds,
                                           Gedim::IMeshDAO &mesh,
                                           std::vector<std::vector<unsigned int>> &meshCell3DsOriginalCell3Ds);

    static void CreateRandomlyDeformedQuadrilaterals(const Gedim::GeometryUtilities &geometryUtilities,
                                                     const Eigen::Vector3d &rectangleOrigin,
                                                     const Eigen::Vector3d &rectangleBaseTangent,
                                                     const Eigen::Vector3d &rectangleHeightTangent,
                                                     const unsigned int &numQuadrilateralsBaseTangent,
                                                     const unsigned int &numQuadrilateralsHeightTangent,
                                                     const double &maxDeformingPercentageBase,
                                                     const double &maxDeformingPercentageHeight,
                                                     Gedim::IMeshDAO &mesh);

    static void CreateDistortedQuadrilaterals(const Gedim::GeometryUtilities &geometryUtilities,
                                              const Eigen::Vector3d &rectangleOrigin,
                                              const Eigen::Vector3d &rectangleBaseTangent,
                                              const Eigen::Vector3d &rectangleHeightTangent,
                                              const unsigned int &numQuadrilateralsBaseTangent,
                                              const unsigned int &numQuadrilateralsHeightTangent,
                                              Gedim::IMeshDAO &mesh);
    /// \brief Given a set of Cell2Ds find the common Cell0Ds
    /// \param cell2DsIndex the cell2Ds index
    /// \param mesh the mesh
    /// \return the Cell0D indices
    static std::vector<unsigned int> FindCell2DsCommonVertices(const std::vector<unsigned int> &cell2DsIndex,
                                                               const Gedim::IMeshDAO &mesh);

    /// \brief Given a set of Cell2Ds find the common Cell1Ds
    /// \param cell2DsIndex the cell2Ds index
    /// \param mesh the mesh
    /// \return the Cell1D indices
    static std::vector<unsigned int> FindCell2DsCommonEdges(const std::vector<unsigned int> &cell2DsIndex,
                                                            const Gedim::IMeshDAO &mesh);

    static Gedim::MeshUtilities::FindConcaveCell3DFacesConvexCell2DResult FindConcaveCell3DFacesConvexCell2D(
        const Gedim::GeometryUtilities &geometryUtilities,
        const unsigned int &concaveCell3DIndex,
        const Gedim::IMeshDAO &mesh,
        const Gedim::IMeshDAO &convexMesh,
        const std::vector<unsigned int> &convexCell3DIndices,
        const std::vector<Eigen::MatrixXd> &concaveCell3DFaces3DVertices,
        const std::vector<Eigen::MatrixXd> &concaveCell3DFaces2DVertices,
        const std::vector<Eigen::Vector3d> &concaveCell3DFacesTranslation,
        const std::vector<Eigen::Matrix3d> &concaveCell3DFacesRotationMatrix,
        const std::vector<Eigen::Vector3d> &concaveCell3DFacesNormal,
        const std::vector<std::vector<Eigen::MatrixXd>> &convexCell3DsFaces3DVertices,
        const std::vector<std::vector<std::vector<unsigned int>>> &convexCell3DsFacesUnalignedVertices);

    static Gedim::MeshUtilities::FindConcaveCell3DFacesConvexCell2DResult FindConcaveCell3DFacesConvexCell2D(
        const Gedim::GeometryUtilities &geometryUtilities,
        const unsigned int &concaveCell3DIndex,
        const Gedim::IMeshDAO &mesh,
        const std::vector<Eigen::MatrixXd> &concaveCell3DTetra,
        const std::vector<std::vector<Eigen::Matrix3d>> &concaveCell3D_faces_2D_triangles,
        const std::vector<Eigen::MatrixXd> &concaveCell3DFaces3DVertices,
        const std::vector<Eigen::Vector3d> &concaveCell3DFacesTranslation,
        const std::vector<Eigen::Matrix3d> &concaveCell3DFacesRotationMatrix,
        const std::vector<Eigen::Vector3d> &concaveCell3DFacesNormal,
        const std::vector<std::vector<Eigen::MatrixXd>> &convexCell3DsFaces3DVertices,
        const std::vector<std::vector<std::vector<unsigned int>>> &convexCell3DsFacesUnalignedVertices);

    static Gedim::MeshUtilities::FindPointMeshPositionResult FindPointMeshPosition(const MeshUtilities::FindPointCell2DResult &find_cell2D_result,
                                                                                   const Gedim::IMeshDAO &mesh);

    static Gedim::MeshUtilities::FindPointMeshPositionResult FindPointMeshPosition(const MeshUtilities::FindPointCell3DResult &find_cell3D_result,
                                                                                   const Gedim::IMeshDAO &mesh);

    static Gedim::MeshUtilities::FindPointCell2DResult FindPointCell2D(const Gedim::GeometryUtilities &geometryUtilities,
                                                                       const Eigen::Vector3d &point,
                                                                       const Gedim::IMeshDAO &mesh,
                                                                       const std::vector<Eigen::MatrixXd> &cell2DsVertices,
                                                                       const std::vector<Eigen::MatrixXd> &cell2DsBoundingBox,
                                                                       const bool find_only_first_cell2D = true,
                                                                       const unsigned int starting_cell2D_index = 0);

    static Gedim::MeshUtilities::FindPointCell3DResult FindPointCell3D(
        const Gedim::GeometryUtilities &geometryUtilities,
        const Eigen::Vector3d &point,
        const Gedim::IMeshDAO &mesh,
        const std::vector<std::vector<Eigen::MatrixXi>> &cell3DsFaces,
        const std::vector<std::vector<Eigen::MatrixXd>> &cell3DsFaceVertices,
        const std::vector<std::vector<Eigen::MatrixXd>> &cell3DsFaceRotatedVertices,
        const std::vector<std::vector<Eigen::Vector3d>> &cell3DsFaceNormals,
        const std::vector<std::vector<bool>> &cell3DsFaceNormalDirections,
        const std::vector<std::vector<Eigen::Vector3d>> &cell3DsFaceTranslations,
        const std::vector<std::vector<Eigen::Matrix3d>> &cell3DsFaceRotationMatrices,
        const std::vector<Eigen::MatrixXd> &cell3DsBoundingBox,
        const bool find_only_first_cell3D = true,
        const unsigned int starting_cell3D_index = 0);

    static Gedim::MeshUtilities::FindPointCell3DResult FindPointCell3D(
        const Gedim::GeometryUtilities &geometryUtilities,
        const Eigen::Vector3d &point,
        const Gedim::IMeshDAO &mesh,
        const std::vector<std::vector<Eigen::MatrixXi>> &cell3DsFaces,
        const std::vector<std::vector<Eigen::MatrixXd>> &cell3DsFaceVertices,
        const std::vector<std::vector<Eigen::MatrixXd>> &cell3DsFaceRotatedVertices,
        const std::vector<std::vector<Eigen::Vector3d>> &cell3DsFaceNormals,
        const std::vector<std::vector<bool>> &cell3DsFaceNormalDirections,
        const std::vector<std::vector<Eigen::Vector3d>> &cell3DsFaceTranslations,
        const std::vector<std::vector<Eigen::Matrix3d>> &cell3DsFaceRotationMatrices,
        const std::vector<Eigen::MatrixXd> &cell3DsBoundingBox,
        const std::vector<std::vector<Eigen::MatrixXd>> &cell3DsTetrahedra,
        const bool find_only_first_cell3D,
        const unsigned int starting_cell3D_index);

    /// \brief Agglomerate Triangles with one vertex in common
    /// \param trianglesIndexToAgglomerate the cell2Ds triangular index in the mesh
    /// \param triangularMesh the triangular mesh
    /// \return the agglomearted polygon indices
    /// \note the triangular index shall be done counterclockwise
    static Gedim::MeshUtilities::AgglomerateTrianglesResult AgglomerateTriangles(const std::vector<unsigned int> &trianglesIndexToAgglomerate,
                                                                                 Gedim::IMeshDAO &triangularMesh);

    static Gedim::MeshUtilities::AgglomerateMeshFromTriangularMeshResult AgglomerateMeshFromTriangularMesh(
        const std::vector<std::vector<unsigned int>> &trianglesIndicesToAgglomerate,
        Gedim::IMeshDAO &triangularMesh);

    /// \brief Import Agglomeration mesh Information From file Csv
    /// \param Gedim::GeometryUtilities the geometry utilities
    /// \param originalMesh the original mesh
    /// \param agglomeratedMesh the agglomerated mesh
    /// \param fileName the csv file name
    /// \param separator the csv file separator
    /// \param originalCell0DToAgglomeratedCell0Ds original Cell0Ds to agglomerated Cell0Ds
    /// \param originalCell1DToAgglomeratedCell1Ds original Cell1Ds to agglomerated Cell1Ds
    /// \param originalCell2DToAgglomeratedCell2Ds original Cell2Ds to agglomerated Cell2Ds
    /// \param agglomeratedCell0DToOriginalCell0Ds agglomerated Cell0Ds to original Cell0Ds
    /// \param agglomeratedCell1DToOriginalCell1Ds agglomerated Cell1Ds to original Cell1Ds
    /// \param agglomeratedCell2DToOriginalCell2Ds agglomerated Cell2Ds to original Cell2Ds
    static Gedim::MeshUtilities::AgglomerationInformation ImportAgglomerationInformationFromCsv(const Gedim::GeometryUtilities &geometryUtilities,
                                                                                                const Gedim::IMeshDAO &originalMesh,
                                                                                                const Gedim::IMeshDAO &agglomeratedMesh,
                                                                                                const std::string &fileName,
                                                                                                const char &separator);

    /// \brief Import Agglomeration mesh Information From file OFF
    /// \param Gedim::GeometryUtilities the geometry utilities
    /// \param originalMesh the original mesh
    /// \param agglomeratedMesh the agglomerated mesh
    /// \param fileName the csv file name
    /// \param separator the csv file separator
    /// \param originalCell0DToAgglomeratedCell0Ds original Cell0Ds to agglomerated Cell0Ds
    /// \param originalCell1DToAgglomeratedCell1Ds original Cell1Ds to agglomerated Cell1Ds
    /// \param originalCell2DToAgglomeratedCell2Ds original Cell2Ds to agglomerated Cell2Ds
    /// \param agglomeratedCell0DToOriginalCell0Ds agglomerated Cell0Ds to original Cell0Ds
    /// \param agglomeratedCell1DToOriginalCell1Ds agglomerated Cell1Ds to original Cell1Ds
    /// \param agglomeratedCell2DToOriginalCell2Ds agglomerated Cell2Ds to original Cell2Ds
    static Gedim::MeshUtilities::AgglomerationInformation ImportAgglomerationInformationFromOFF(const Gedim::GeometryUtilities &geometryUtilities,
                                                                                                const Gedim::IMeshDAO &originalMesh,
                                                                                                const Gedim::IMeshDAO &agglomeratedMesh,
                                                                                                const std::string &fileName,
                                                                                                const char &separator);

    /// \brief Export mesh to csv file
    static void ExportMeshToCsv(const Gedim::IMeshDAO &mesh, const char &separator, const std::string &exportFolderPath);

    /// \brief Export 2D concave mesh to csv file
    static void ExportConcaveMesh2DToCsv(const Gedim::IMeshDAO &mesh,
                                         const std::vector<std::vector<unsigned int>> &convexCell2DsIndex,
                                         const char &separator,
                                         const std::string &exportFolderPath);

    static std::vector<int> MarkCells(const std::function<Eigen::VectorXi(const Eigen::MatrixXd &)> &marking_function,
                                      const std::vector<Eigen::MatrixXd> &cells_points,
                                      const unsigned int default_mark);

    static Gedim::MeshUtilities::Intersect_mesh_polyhedron_result Intersect_mesh_polyhedron(
        const Gedim::GeometryUtilities &geometry_utilities,
        const Eigen::MatrixXd &polyhedron_vertices,
        const Eigen::MatrixXi &polyhedron_edges,
        const std::vector<Eigen::MatrixXd> &polyhedron_edges_vertices,
        const Eigen::MatrixXd &polyhedron_edges_tangent,
        const std::vector<Eigen::MatrixXd> &polyhedron_edges_boudingBox,
        const std::vector<Eigen::MatrixXi> &polyhedron_faces,
        const std::vector<Eigen::MatrixXd> &polyhedron_faces_vertices,
        const std::vector<Eigen::MatrixXd> &polyhedron_faces_rotated_vertices,
        const std::vector<Eigen::Vector3d> &polyhedron_faces_normals,
        const std::vector<bool> &polyhedron_faces_normal_direction,
        const std::vector<Eigen::Vector3d> &polyhedron_faces_translation,
        const std::vector<Eigen::Matrix3d> &polyhedron_faces_rotation_matrix,
        const std::vector<Eigen::MatrixXd> &polyhedron_faces_boudingBox,
        const Eigen::MatrixXd &polyhedron_boudingBox,
        const Gedim::IMeshDAO &mesh,
        const std::vector<Eigen::MatrixXd> &mesh_cell1Ds_boudingBox,
        const std::vector<Eigen::MatrixXd> &mesh_cell1Ds_vertices,
        const std::vector<Eigen::Vector3d> &mesh_cell1Ds_tangent,
        const std::vector<Eigen::MatrixXd> &mesh_cell2Ds_vertices,
        const std::vector<Eigen::Vector3d> &mesh_cell2Ds_normal,
        const std::vector<Eigen::MatrixXd> &mesh_cell2Ds_2D_vertices,
        const std::vector<Eigen::Vector3d> &mesh_cell2Ds_translation,
        const std::vector<Eigen::Matrix3d> &mesh_cell2Ds_rotation_matrix,
        const std::vector<Eigen::MatrixXd> &mesh_cell2Ds_boudingBox,
        const std::vector<Eigen::MatrixXd> &mesh_cell3Ds_boudingBox,
        const std::vector<std::vector<Eigen::MatrixXi>> &mesh_cell3Ds_faces,
        const std::vector<std::vector<Eigen::MatrixXd>> &mesh_cell3Ds_faces_vertices,
        const std::vector<std::vector<Eigen::MatrixXd>> &mesh_cell3Ds_faces_2D_vertices,
        const std::vector<std::vector<Eigen::Vector3d>> &mesh_cell3Ds_faces_normal,
        const std::vector<std::vector<bool>> &mesh_cell3Ds_faces_normal_directions,
        const std::vector<std::vector<Eigen::Vector3d>> &mesh_cell3Ds_faces_translation,
        const std::vector<std::vector<Eigen::Matrix3d>> &mesh_cell3Ds_faces_rotation_matrix);

    static void SetPolygonMeshMarkers(const Gedim::GeometryUtilities &geometryUtilities,
                                      const Eigen::MatrixXd &polygonVertices,
                                      const std::vector<unsigned int> &cell0DMarkers,
                                      const std::vector<unsigned int> &cell1DMarkers,
                                      Gedim::IMeshDAO &mesh);

    static void CheckMeshGeometricData2D(const Gedim::MeshUtilities::CheckMeshGeometricData2DConfiguration &configuration,
                                         const Gedim::GeometryUtilities &geometryUtilities,
                                         const Gedim::IMeshDAO &mesh,
                                         const MeshGeometricData2D &geometricData);

    static void CreatePolygonIntersectionMesh(const Gedim::GeometryUtilities &geometry_utilities,
                                              const Eigen::MatrixXd &interface_vertices,
                                              Gedim::IMeshDAO &mesh);

    static bool CollapseCell1D(const unsigned int cell1D_index, Gedim::IMeshDAO &mesh);

    static void Mesh3DFromPolyhedra(const Gedim::GeometryUtilities &geometry_utilities,
                                    const std::vector<Gedim::GeometryUtilities::Polyhedron> &polyhedra,
                                    Gedim::IMeshDAO &mesh);

    static MarchingCubesData MarchingCubes(const Gedim::GeometryUtilities &geometry_utilities,
                                           const std::function<Eigen::VectorXd(const Eigen::MatrixXd &)> &surface_implicit_equation,
                                           const Gedim::IMeshDAO &mesh);
};

} // namespace Gedim

#endif // __MeshUtilities_H
