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

#include "CommonUtilities.hpp"
#include "MeshUtilities.hpp"

namespace Gedim
{
// ***************************************************************************
using Segment = std::array<unsigned int, 2>;

inline constexpr unsigned int None = std::numeric_limits<unsigned int>::max();

inline constexpr std::array<std::array<int, 4>, 16> SegmentTable = {{
    {{-1, -1, -1, -1}}, // 0
    {{3, 0, -1, -1}},   // 1
    {{0, 1, -1, -1}},   // 2
    {{3, 1, -1, -1}},   // 3
    {{1, 2, -1, -1}},   // 4
    {{3, 0, 1, 2}},     // 5: ambiguous
    {{0, 2, -1, -1}},   // 6
    {{3, 2, -1, -1}},   // 7
    {{2, 3, -1, -1}},   // 8
    {{0, 2, -1, -1}},   // 9
    {{0, 1, 2, 3}},     // 10: ambiguous
    {{1, 2, -1, -1}},   // 11
    {{1, 3, -1, -1}},   // 12
    {{0, 1, -1, -1}},   // 13
    {{0, 3, -1, -1}},   // 14
    {{-1, -1, -1, -1}}  // 15
}};

inline constexpr std::array<std::array<int, 2>, 4> SquareEdgeEnds = {{{{0, 1}}, {{1, 2}}, {{2, 3}}, {{3, 0}}}};

inline constexpr std::array<std::array<int, 2>, 4> SquareBits = {{{{0, 0}}, {{1, 0}}, {{1, 1}}, {{0, 1}}}};

struct Square
{
    std::array<unsigned int, 4> vertices;
    std::array<unsigned int, 4> edges;
    Eigen::Vector3d lo, hi;
};
// ***************************************************************************
inline Square BuildSquare(const Gedim::GeometryUtilities &geometry_utilities, const Gedim::IMeshDAO &mesh, unsigned int cell)
{
    if (mesh.Cell2DNumberVertices(cell) != 4 || mesh.Cell2DNumberEdges(cell) != 4)
        throw std::runtime_error("Expected quadrilateral with 4 vertices and 4 edges");

    Square s;
    s.vertices.fill(None);
    s.edges.fill(None);

    s.lo.setConstant(std::numeric_limits<double>::infinity());
    s.hi.setConstant(-std::numeric_limits<double>::infinity());

    for (unsigned int k = 0; k < 4; ++k)
    {
        const auto coordinate = mesh.Cell0DCoordinates(mesh.Cell2DVertex(cell, k));

        for (unsigned int a = 0; a < 3; ++a)
        {
            s.lo[a] = std::min(s.lo[a], coordinate[a]);
            s.hi[a] = std::max(s.hi[a], coordinate[a]);
        }
    }

    const Eigen::Vector3d extent = s.hi - s.lo;

    if (geometry_utilities.IsValueLowerOrEqual(extent.x(), 0.0, geometry_utilities.Tolerance1D()) ||
        geometry_utilities.IsValueLowerOrEqual(extent.y(), 0.0, geometry_utilities.Tolerance1D()) ||
        !geometry_utilities.IsValueZero(extent.z(), geometry_utilities.Tolerance1D()))
        throw std::runtime_error("Expected nondegenerate axis-aligned XY quadrilateral");

    for (unsigned int k = 0; k < 4; ++k)
    {
        const unsigned int vertex = mesh.Cell2DVertex(cell, k);

        const auto coordinate = mesh.Cell0DCoordinates(vertex);

        int bits[2];

        for (unsigned int a = 0; a < 2; ++a)
        {
            const double t = (coordinate[a] - s.lo[a]) / extent[a];

            if (geometry_utilities.IsValueZero(t, geometry_utilities.Tolerance1D()))
                bits[a] = 0;
            else if (geometry_utilities.IsValueZero(t - 1.0, geometry_utilities.Tolerance1D()))
                bits[a] = 1;
            else
                throw std::runtime_error("Cell is not axis-aligned Cartesian");
        }

        int local = -1;

        for (unsigned int j = 0; j < 4; ++j)
        {
            if (SquareBits[j][0] == bits[0] && SquareBits[j][1] == bits[1])
            {
                local = static_cast<int>(j);
                break;
            }
        }

        if (local < 0 || s.vertices[local] != None)
            throw std::runtime_error("Invalid square vertex ordering");

        s.vertices[local] = vertex;
    }

    for (unsigned int e = 0; e < 4; ++e)
    {
        const unsigned int a = s.vertices[SquareEdgeEnds[e][0]];
        const unsigned int b = s.vertices[SquareEdgeEnds[e][1]];

        for (unsigned int j = 0; j < 4; ++j)
        {
            const unsigned int edge = mesh.Cell2DEdge(cell, j);

            const unsigned int origin = mesh.Cell1DOrigin(edge);
            const unsigned int end = mesh.Cell1DEnd(edge);

            if ((origin == a && end == b) || (origin == b && end == a))
            {
                s.edges[e] = edge;
                break;
            }
        }

        if (s.edges[e] == None)
            throw std::runtime_error("Missing square edge");
    }

    return s;
}
// ***************************************************************************
inline Eigen::Vector3d GradientAt(const Square &square, const std::array<double, 4> &values, const Eigen::Vector3d &point)
{
    const Eigen::Vector3d extent = square.hi - square.lo;

    const double x = (point.x() - square.lo.x()) / extent.x();
    const double y = (point.y() - square.lo.y()) / extent.y();

    Eigen::Vector3d gradient = Eigen::Vector3d::Zero();

    for (unsigned int i = 0; i < 4; ++i)
    {
        const double wx = (SquareBits[i][0] ? 1.0 : -1.0) / extent.x();

        const double wy = (SquareBits[i][1] ? 1.0 : -1.0) / extent.y();

        gradient.x() += values[i] * wx * (SquareBits[i][1] ? y : (1.0 - y));

        gradient.y() += values[i] * wy * (SquareBits[i][0] ? x : (1.0 - x));
    }

    return gradient;
}
// ***************************************************************************
struct SegmentHash
{
    std::size_t operator()(const Segment &a) const noexcept
    {
        std::size_t h = 0;

        for (auto x : a)
            h ^= std::hash<unsigned int>{}(x) + 0x9e3779b9 + (h << 6) + (h >> 2);

        return h;
    }
};
// ***************************************************************************
inline MeshUtilities::MarchingSquaresResult::SurfaceSegmentsData CreateSegments(
    const Gedim::GeometryUtilities &geometry_utilities,
    const Eigen::VectorXd &cell0Ds_surface_value,
    const std::vector<Gedim::MeshUtilities::PointSurfacePositions> &cell0Ds_surface_position,
    const std::vector<Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::SurfacePointData> &surface_intersection_points,
    const std::vector<unsigned int> &cell0Ds_surface_points,
    const std::vector<unsigned int> &cell1Ds_surface_points,
    const Gedim::IMeshDAO &mesh)
{
    using Position = Gedim::MeshUtilities::PointSurfacePositions;

    if (static_cast<unsigned int>(cell0Ds_surface_value.size()) != mesh.Cell0DTotalNumber() ||
        cell0Ds_surface_position.size() != mesh.Cell0DTotalNumber() ||
        cell0Ds_surface_points.size() != mesh.Cell0DTotalNumber() || cell1Ds_surface_points.size() != mesh.Cell1DTotalNumber())
        throw std::invalid_argument("Surface arrays do not match mesh dimensions");

    MeshUtilities::MarchingSquaresResult::SurfaceSegmentsData result;

    std::unordered_set<Segment, SegmentHash> seen;

    result.cell2Ds_surface_segments.resize(mesh.Cell2DTotalNumber());

    result.surface_segments.reserve(2 * mesh.Cell2DTotalNumber());

    for (unsigned int cell = 0; cell < mesh.Cell2DTotalNumber(); ++cell)
    {
        result.cell2Ds_surface_segments.at(cell).reserve(2);

        const Square s = BuildSquare(geometry_utilities, mesh, cell);

        unsigned int mask = 0;
        std::array<double, 4> values;

        for (unsigned int i = 0; i < 4; ++i)
        {
            const unsigned int v = s.vertices[i];

            values[i] = cell0Ds_surface_value[v];

            if (cell0Ds_surface_position[v] == Position::Negative)
                mask |= 1u << i;
        }

        // OnSurface belongs to the nonnegative side.
        if (mask == 0 || mask == 15)
            continue;

        auto row = SegmentTable[mask];

        // Resolve ambiguous saddle cases using the
        // bilinear asymptotic decider.
        if (mask == 5 || mask == 10)
        {
            const double q = values[0] * values[2] - values[1] * values[3];

            // q > 0: connect (0,1) and (2,3).
            // q < 0: connect (3,0) and (1,2).
            //
            // q == 0: deterministic tie-breaking.
            if (q > 0.0)
                row = {{0, 1, 2, 3}};
            else if (q < 0.0)
                row = {{3, 0, 1, 2}};
            else
                row = SegmentTable[mask];
        }

        for (unsigned int j = 0; j < 4 && row[j] >= 0; j += 2)
        {
            if (j + 1 >= 4 || row[j + 1] < 0)
                throw std::runtime_error("Malformed segment table");

            Segment segment{};

            for (unsigned int k = 0; k < 2; ++k)
            {
                const unsigned int e = static_cast<unsigned int>(row[j + k]);

                const unsigned int a = s.vertices[SquareEdgeEnds[e][0]];
                const unsigned int b = s.vertices[SquareEdgeEnds[e][1]];

                const bool az = cell0Ds_surface_position[a] == Position::OnSurface;

                const bool bz = cell0Ds_surface_position[b] == Position::OnSurface;

                // Same convention as Marching Cubes.
                const unsigned int id = (az && bz) ? cell0Ds_surface_points[std::min(a, b)]
                                        : az       ? cell0Ds_surface_points[a]
                                        : bz       ? cell0Ds_surface_points[b]
                                                   : cell1Ds_surface_points[s.edges[e]];

                if (id == None || id >= surface_intersection_points.size())
                    throw std::runtime_error("Required precomputed surface point is missing");

                segment[k] = id;
            }

            // Remove zero-length segments by topology.
            if (segment[0] == segment[1])
                continue;

            const auto &p0 = surface_intersection_points[segment[0]].coordinate;

            const auto &p1 = surface_intersection_points[segment[1]].coordinate;

            if (p0.size() < 3 || p1.size() < 3)
                throw std::runtime_error("Expected 3D intersection points");

            const Eigen::Vector3d a = p0.head<3>();
            const Eigen::Vector3d b = p1.head<3>();

            const Eigen::Vector3d tangent = b - a;

            const double length2 = tangent.squaredNorm();

            if (geometry_utilities.IsValueZero(length2, geometry_utilities.Tolerance1DSquared()))
                continue;

            // Orient segment consistently with the
            // bilinear scalar field gradient.
            const Eigen::Vector3d gradient = GradientAt(s, values, (a + b) / 2.0);

            // For tangent (tx,ty), its right normal is
            // (ty,-tx). Align it with the positive gradient.
            const double orientation = tangent.y() * gradient.x() - tangent.x() * gradient.y();

            if (geometry_utilities.IsValuePositive(gradient.squaredNorm(), geometry_utilities.Tolerance1DSquared()) &&
                geometry_utilities.IsValueNegative(orientation, geometry_utilities.Tolerance1D()))
                std::swap(segment[0], segment[1]);

            Segment key = segment;
            std::sort(key.begin(), key.end());

            if (seen.insert(key).second)
            {
                result.cell2Ds_surface_segments.at(cell).push_back(result.surface_segments.size());

                result.surface_segments.push_back({segment, cell});
            }
        }
    }

    for (auto &cell_segments : result.cell2Ds_surface_segments)
        cell_segments.shrink_to_fit();

    result.cell2Ds_surface_segments.shrink_to_fit();
    result.surface_segments.shrink_to_fit();

    return result;
}
// ***************************************************************************
Gedim::MeshUtilities::MarchingSquaresResult Gedim::MeshUtilities::MarchingSquares(
    const Gedim::GeometryUtilities &geometry_utilities,
    const std::function<Eigen::VectorXd(const Eigen::MatrixXd &)> &surface_implicit_equation,
    const Gedim::IMeshDAO &mesh)
{
    using SurfacePosition = Gedim::MeshUtilities::PointSurfacePositions;

    Gedim::MeshUtilities::MarchingSquaresResult result;

    result.cell0Ds_surface_value = surface_implicit_equation(mesh.Cell0DsCoordinates());

    {
        auto marking_function = [&geometry_utilities](const Eigen::VectorXd &points_surface_value) {
            std::vector<SurfacePosition> marks(points_surface_value.size());

            for (unsigned int p = 0; p < points_surface_value.size(); p++)
            {
                const auto surface_value = points_surface_value[p];

                if (geometry_utilities.IsValuePositive(surface_value, geometry_utilities.Tolerance1D()))
                {
                    marks[p] = SurfacePosition::Positive;
                }
                else if (geometry_utilities.IsValueZero(surface_value, geometry_utilities.Tolerance1D()))
                    marks[p] = SurfacePosition::OnSurface;
                else if (geometry_utilities.IsValueNegative(surface_value, geometry_utilities.Tolerance1D()))
                    marks[p] = SurfacePosition::Negative;
                else
                    throw std::runtime_error("unknown cell0D position");
            }

            return marks;
        };

        result.cell0Ds_surface_position = marking_function(result.cell0Ds_surface_value);
    }

    result.surface_intersections =
        IntersectMeshWithSurface(geometry_utilities, result.cell0Ds_surface_value, result.cell0Ds_surface_position, mesh);

    result.surface_segments = CreateSegments(geometry_utilities,
                                             result.cell0Ds_surface_value,
                                             result.cell0Ds_surface_position,
                                             result.surface_intersections.surface_intersection_points,
                                             result.surface_intersections.cell0Ds_surface_points,
                                             result.surface_intersections.cell1Ds_surface_points,
                                             mesh);

    return result;
}
// ***************************************************************************
} // namespace Gedim
