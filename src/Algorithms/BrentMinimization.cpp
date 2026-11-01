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

#include "BrentMinimization.hpp"

#include <iostream>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <algorithm>
#include <numbers>


namespace Gedim
{
  namespace Algorithms
  {
    namespace Minimization
    {
      constexpr double golden_section = 2.0 - std::numbers::phi_v<double>;
      // ***************************************************************************
      BrentResult brent_minimize(
          const std::function<double(double)>& f,
          const double& a,
          const double& b,
          const double& tol,
          const unsigned int max_iter)
      {
        if (!(a < b) || !std::isfinite(a) || !std::isfinite(b) ||
            !(tol > 0.0) || !std::isfinite(tol) || max_iter <= 0) {
          throw std::invalid_argument("Invalid Brent parameters");
        }

        // Initial point
        double x = a + golden_section * (b - a);
        double w = x;
        double v = x;

        double fx = f(x);
        double fw = fx;
        double fv = fx;

        // Step history
        double d = 0.0;
        double e = 0.0;

        double left = a;
        double right = b;

        for (unsigned int iter = 1; iter <= max_iter; ++iter)
        {
          const double midpoint = 0.5 * (right + left);

          const double abs_tol = std::abs(x) < 1.0 ? tol :
                                                     tol * std::abs(x);
          if (std::abs(x - midpoint) <= abs_tol &&
              right <= 2.0 * abs_tol + left) // max(|x-m|, 0.5 (r-l)) / x < tol
          {
            return { x, fx, iter - 1, true };
          }

          bool use_golden = true;

          // Attempt parabolic interpolation
          if (std::abs(e) > abs_tol) {

            double r = (x - w) * (fx - fv);
            double q = (x - v) * (fx - fw);
            double p = (x - v) * q - (x - w) * r;

            q = 2.0 * (q - r);

            if (q > 0.0)
              p = -p;

            q = std::abs(q);

            const double previous_e = e;
            e = d;

            // Accept parabolic step only if safe
            if (q > 0.0 &&
                std::abs(p) < std::abs(0.5 * q * previous_e) &&
                p > q * (left - x) &&
                p < q * (right - x)) {

              d = p / q;
              const double u = x + d;

              // Avoid points too close to boundaries
              if ((u - left) < abs_tol || (right - u) < abs_tol) {
                d = std::copysign(abs_tol, midpoint - x);
              }

              use_golden = false;
            }
          }

          // Golden-section fallback
          if (use_golden)
          {
            e = (x < midpoint) ? (right - x) : (left - x);
            d = golden_section * e;
          }

          const double u = x + (
                             std::abs(d) >= abs_tol
                             ? d
                             : std::copysign(abs_tol, d)
                               );

          const double fu = f(u);

          // Update interval and interpolation points
          if (fu <= fx)
          {
            if (u < x)
              right = x;
            else
              left = x;

            v = w;
            fv = fw;

            w = x;
            fw = fx;

            x = u;
            fx = fu;

          }
          else
          {

            if (u < x)
              left = u;
            else
              right = u;

            if (fu <= fw || w == x)
            {
              v = w;
              fv = fw;

              w = u;
              fw = fu;

            }
            else if (fu <= fv || v == x || v == w)
            {
              v = u;
              fv = fu;
            }
          }
        }

        return { x, fx, max_iter, false };
      }
      // ***************************************************************************
    } // namespace Minimization
  } // namespace Algorithms
} // namespace Gedim
