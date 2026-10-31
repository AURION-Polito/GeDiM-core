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


namespace Gedim
{
  namespace Algorithms
  {
    namespace Minimization
    {
      constexpr double golden_section = 0.3819660112501051;
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

        const double eps = std::numeric_limits<double>::epsilon();

        double left = a;
        double right = b;

        for (unsigned int iter = 1; iter <= max_iter; ++iter) {

          const double midpoint = left + 0.5 * (right - left);

          const double tol1 = tol * std::abs(x) + eps;
          const double tol2 = 2.0 * tol1;

          // Convergence test
          if (std::abs(x - midpoint) <=
              tol2 - 0.5 * (right - left)) {
            return {x, fx, iter - 1, true};
          }

          bool use_golden = true;

          // Attempt parabolic interpolation
          if (std::abs(e) > tol1) {

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
              if ((u - left) < tol2 || (right - u) < tol2) {
                d = std::copysign(tol1, midpoint - x);
              }

              use_golden = false;
            }
          }

          // Golden-section fallback
          if (use_golden) {
            e = (x < midpoint) ? (right - x) : (left - x);
            d = golden_section * e;
          }

          // Minimum step size
          const double u = x + (
                             std::abs(d) >= tol1
                             ? d
                             : std::copysign(tol1, d)
                               );

          const double fu = f(u);

          // Update interval and interpolation points
          if (fu <= fx) {

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

          } else {

            if (u < x)
              left = u;
            else
              right = u;

            if (fu <= fw || w == x) {
              v = w;
              fv = fw;

              w = u;
              fw = fu;

            } else if (fu <= fv || v == x || v == w) {
              v = u;
              fv = fu;
            }
          }
        }

        return {x, fx, max_iter, false};
      }
      // ***************************************************************************
    } // namespace Minimization
  } // namespace Algorithms
} // namespace Gedim
