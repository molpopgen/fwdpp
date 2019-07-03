#ifndef FWDPP_POISSON_POINT_HPP__
#define FWDPP_POISSON_POINT_HPP__

#include <cmath>
#include <stdexcept>
#include <gsl/gsl_randist.h>
#include "genetic_map_unit.hpp"

namespace fwdpp
{
    struct poisson_point : public genetic_map_unit
    {
        const double position, mean;
        poisson_point(const double pos, const double m)
            : genetic_map_unit(), position(pos), mean(m)
        {
            if (!std::isfinite(pos))
                {
                    throw std::invalid_argument("position must be finite");
                }
            if (!std::isfinite(mean))
                {
                    throw std::invalid_argument("mean must be finite");
                }
            if (mean < 0)
                {
                    throw std::invalid_argument("mean must be non-negative");
                }
        }

        void
        operator()(const gsl_rng* r,
                   std::vector<double>& breakpoints) const final
        {
            unsigned n = gsl_ran_poisson(r, mean);
            if (n % 2 != 0.0)
                {
                    breakpoints.push_back(position);
                }
        }

        std::unique_ptr<genetic_map_unit>
        clone() const final
        {
            return std::unique_ptr<genetic_map_unit>(
                new poisson_point(position, mean));
        }
    };
} // namespace fwdpp

#endif
