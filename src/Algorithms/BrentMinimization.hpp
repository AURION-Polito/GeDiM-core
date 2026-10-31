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

#ifndef __Brent_Minimization_H
#define __Brent_Minimization_H

#include <functional>

namespace Gedim
{
  namespace Algorithms
  {
    namespace Minimization
    {
      struct BrentResult final
      {
          double min;
          double f_min;
          unsigned int iterations;
          bool converged;
      };

      BrentResult brent_minimize(const std::function<double(double)>& f,
                                 const double& a,
                                 const double& b,
                                 const double& tol,
                                 const unsigned int max_iter);
    } // namespace Minimization
  } // namespace Algorithms
} // namespace Gedim

#endif // __Fast_Marching_Method_H
