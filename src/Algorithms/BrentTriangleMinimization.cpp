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
#include "BrentSegmentMinimization.hpp"

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

        const double max_length = std::max({ std::abs(a.at(1) - a.at(0)),
                                             std::abs(b.at(1) - b.at(0)),
                                             std::abs(c.at(1) - c.at(0)) });

        const double abs_tol = std::abs(max_length) < 1.0 ? tol :
                                                            tol * max_length;

        double s_min = 0.0;

        auto g =
            [&](double t)
        {
          if (std::abs(t) < abs_tol)
            return f(a.at(0), a.at(1));

          auto h =
              [&](double s)
          {
            const double P_0 = (1.0 - t) * a.at(0) + t * ((1.0 - s) * b.at(0) + s * c.at(0));
            const double P_1 = (1.0 - t) * a.at(1) + t * ((1.0 - s) * b.at(1) + s * c.at(1));

            return f(P_0, P_1);
          };

          const auto s_min_res = brent_segment_minimization(h, 0.0, 1.0, tol, max_iter);

          if (!s_min_res.converged)
            return std::numeric_limits<double>::quiet_NaN();

          s_min = s_min_res.min.at(0);
          const double P_s_0 = (1.0 - t) * a.at(0) + t * ((1.0 - s_min) * b.at(0) + s_min * c.at(0));
          const double P_s_1 = (1.0 - t) * a.at(1) + t * ((1.0 - s_min) * b.at(1) + s_min * c.at(1));

          return f(P_s_0, P_s_1);
        };

        const auto t_min_res = brent_segment_minimization(g, 0.0, 1.0, tol, max_iter);
        if (std::isnan(t_min_res.f_min))
          return { { std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN() }, std::numeric_limits<double>::quiet_NaN(), max_iter, false };

        double t_min = t_min_res.min.at(0);

        const double P_st_0 = (1.0 - t_min) * a.at(0) + t_min * ((1.0 - s_min) * b.at(0) + s_min * c.at(0));
        const double P_st_1 = (1.0 - t_min) * a.at(1) + t_min * ((1.0 - s_min) * b.at(1) + s_min * c.at(1));

        return { { P_st_0, P_st_1 }, t_min_res.f_min, t_min_res.iterations, t_min_res.converged };
      }
      // ***************************************************************************
    } // namespace Minimization
  } // namespace Algorithms
} // namespace Gedim
