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
        detail::genetic_map_unit_cast_function cast;
        const bool discrete_;

        [[deprecated("this constructor is deprecated as of 0.9.0")]] binomial_point(
            const double pos, const double pr)
            : position(pos), prob(pr),
              cast(detail::generate_genetic_map_unit_cast_function(false)),
              discrete_(false)

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
                    throw std::invalid_argument("probability must be 0 <= x <= 1");
                }
        }

        binomial_point(const double pos, const double pr, bool discrete)
            : position(pos), prob(pr),
              cast(detail::generate_genetic_map_unit_cast_function(discrete)),
              discrete_(discrete)
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
                    throw std::invalid_argument("probability must be 0 <= x <= 1");
                }
        }
        void
        operator()(const gsl_rng* r, std::vector<double>& breakpoints) const final
        {
            if (gsl_rng_uniform(r) <= prob)
                {
                    breakpoints.push_back(cast(position));
                }
        }

        std::unique_ptr<genetic_map_unit>
        clone() const final
        {
            return std::unique_ptr<genetic_map_unit>(
                new binomial_point(position, prob, discrete_));
        }

        virtual bool
        discrete() const final
        {
            return discrete_;
        }
    };
} // namespace fwdpp

#endif

