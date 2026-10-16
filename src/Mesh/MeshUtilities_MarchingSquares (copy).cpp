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
  Gedim::MeshUtilities::CreateVerticesAndEdgesFromMeshSurfaceIntersectionResult Gedim::MeshUtilities::CreateVerticesAndEdgesFromMeshSurfaceIntersection(const std::vector<Gedim::MeshUtilities::PointSurfacePositions>& cell0Ds_surface_position,
                                                                                                                                                        const Gedim::MeshUtilities::IntersectMeshWithSurfaceResult& surface_intersections,
                                                                                                                                                        const Gedim::IMeshDAO &mesh)
  {

  }
  // ***************************************************************************
  Gedim::MeshUtilities::ConformMesh2DToSurfaceResult Gedim::MeshUtilities::ConformMesh2DToSurface(
      const Gedim::IMeshDAO &mesh)
  {
    Gedim::MeshUtilities::ConformMesh2DToSurfaceResult result;
    
    result.cell0Ds_cell1Ds = CreateVerticesAndEdgesFromMeshSurfaceIntersection();

    return result;
  }
  // ***************************************************************************
} // namespace Gedim
