#ifndef FWDPP_BINOMIAL_INTERVAL_HPP__
#define FWDPP_BINOMIAL_INTERVAL_HPP__

#include <cmath>
#include <stdexcept>
#include <gsl/gsl_randist.h>
#include "genetic_map_unit.hpp"

namespace fwdpp
{
    struct binomial_interval : public genetic_map_unit
    {
        const double beg, end, prob;
        binomial_interval(double b, double e, double p)
            : genetic_map_unit(), beg(b), end(e), prob(p)
        {
            if (!std::isfinite(beg))
                {
                    throw std::invalid_argument("beg must be finite");
                }
            if (!std::isfinite(end))
                {
                    throw std::invalid_argument("end must be finite");
                }
            if (!std::isfinite(prob))
                {
                    throw std::invalid_argument("prob must be finite");
                }
            if (end <= beg)
                {
                    throw std::invalid_argument(
                        "end must be greater than beg");
                }
            if (prob < 0)
                {
                    throw std::invalid_argument("prob must be non-negative");
                }
        }

        void
        operator()(const gsl_rng* r,
                   std::vector<double>& breakpoints) const final
        {
            if (gsl_rng_uniform(r) <= prob)
                {
                    breakpoints.push_back(gsl_ran_flat(r, beg, end));
                }
        }

        std::unique_ptr<genetic_map_unit>
        clone() const final
        {
            return std::unique_ptr<genetic_map_unit>(
                new binomial_interval(beg, end, prob));
        }
    };
} // namespace fwdpp

#endif
