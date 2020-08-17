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
        detail::genetic_map_unit_cast_function cast;
        const bool discrete_;

        [[deprecated("this constructor is deprecated as of 0.9.0")]] poisson_point(
            const double pos, const double m)
            : genetic_map_unit(), position(pos), mean(m),
              cast(detail::generate_genetic_map_unit_cast_function(false)),
              discrete_(false)
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

        poisson_point(const double pos, const double m, bool discrete)
            : genetic_map_unit(), position(pos), mean(m),
              cast(detail::generate_genetic_map_unit_cast_function(discrete)),
              discrete_(discrete)
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
        operator()(const gsl_rng* r, std::vector<double>& breakpoints) const final
        {
            unsigned n = gsl_ran_poisson(r, mean);
            if (n % 2 != 0.0)
                {
                    breakpoints.push_back(cast(position));
                }
        }

        std::unique_ptr<genetic_map_unit>
        clone() const final
        {
            return std::make_unique<poisson_point>(position, mean, discrete_);
        }

        virtual bool
        discrete() const final
        {
            return discrete_;
        }
    };
} // namespace fwdpp

#endif
