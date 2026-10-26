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

    return result;
  }
  // ***************************************************************************
} // namespace Gedim
