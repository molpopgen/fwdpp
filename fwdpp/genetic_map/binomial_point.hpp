#ifndef FWDPP_BINOMIAL_POINT_HPP
#define FWDPP_BINOMIAL_POINT_HPP

#include <cmath>
#include <stdexcept>
#include <gsl/gsl_randist.h>
#include "genetic_map_unit.hpp"

namespace fwdpp
{
    struct binomial_point : public genetic_map_unit
    {
        double position, prob;
        binomial_point(const double pos, const double pr)
            : position(pos), prob(pr)
        {
            if (!std::isfinite(pos))
                {
                    throw std::invalid_argument("position must be finite");
                }
            if (!std::isfinite(prob))
                {
                    throw std::invalid_argument("probability must be finite");
                }
            if (pr < 0 || pr > 1.)
                {
                    throw std::invalid_argument(
                        "probability must be 0 <= x <= 1");
                }
        }

        void
        operator()(const gsl_rng* r,
                   std::vector<double>& breakpoints) const final
        {
            if (gsl_rng_uniform(r) <= prob)
                {
                    breakpoints.push_back(position);
                }
        }

        std::unique_ptr<genetic_map_unit>
        clone() const final
        {
            return std::unique_ptr<genetic_map_unit>(
                new binomial_point(position, prob));
        }
    };
} // namespace fwdpp

#endif

