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

#ifndef __Minimization_Result_H
#define __Minimization_Result_H

#include <vector>

namespace Gedim
{
  namespace Algorithms
  {
    namespace Minimization
    {
      struct MinimizationResult final
      {
          std::vector<double> min;
          double f_min;
          unsigned int iterations;
          bool converged;
      };
    } // namespace Minimization
  } // namespace Algorithms
} // namespace Gedim

#endif // __Minimization_Result_H
