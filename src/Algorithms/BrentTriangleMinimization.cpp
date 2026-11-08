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

#include "BrentTriangleMinimization.hpp"

#include <cmath>
#include <numbers>
#include <stdexcept>

namespace Gedim
{
  namespace Algorithms
  {
    namespace Minimization
    {
      // ***************************************************************************
      MinimizationResult brent_triangle_minimization(
          const std::function<double(double, double)>& f,
          const std::array<double, 2>& a,
          const std::array<double, 2>& b,
          const std::array<double, 2>& c,
          const double& tol,
          const unsigned int max_iter)
      {
        if (!(tol > 0.0) || !std::isfinite(tol) || max_iter <= 0) {
          throw std::runtime_error("Invalid Brent parameters");
        }

        return { { 0.0, 0.0 }, 0.0, max_iter, false };
      }
      // ***************************************************************************
    } // namespace Minimization
  } // namespace Algorithms
} // namespace Gedim
