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

#ifndef __Brent_Triangle_Minimization_H
#define __Brent_Triangle_Minimization_H

#include "MinimizationResult.hpp"
#include <functional>

namespace Gedim
{
  namespace Algorithms
  {
    namespace Minimization
    {
      MinimizationResult brent_triangle_minimization(const std::function<double(double, double)>& f,
                                 const std::array<double, 2>& a,
                                 const std::array<double, 2>& b,
                                 const std::array<double, 2>& c,
                                 const double& tol,
                                 const unsigned int max_iter);
    } // namespace Minimization
  } // namespace Algorithms
} // namespace Gedim

#endif // __Brent_Triangle_Minimization_H
