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
  Gedim::MeshUtilities::ConformMeshCell2DsToSurfaceData Gedim::MeshUtilities::ConformMeshCell2DsToSurface(
      const Gedim::MeshUtilities::MeshSurfaceIntersectionGraphData& surface_intersection_graph,
      const Gedim::MeshUtilities::MarchingCubesResult::SurfaceTrianglesData& surface_triangles,
      const Gedim::IMeshDAO &mesh)
  {
    Gedim::MeshUtilities::ConformMeshCell2DsToSurfaceData result;

    return result;
  }
  // ***************************************************************************
} // namespace Gedim
