#ifndef FWDPP_POISSON_INTERVAL_HPP__
#define FWDPP_POISSON_INTERVAL_HPP__

#include <cmath>
#include <stdexcept>
#include <gsl/gsl_randist.h>
#include "genetic_map_unit.hpp"

namespace fwdpp
{
    struct poisson_interval : public genetic_map_unit
    {
        const double beg, end, mean;
        detail::genetic_map_unit_cast_function cast;
        const bool discrete_;

        [[deprecated("this constructor is deprecated as of 0.9.0")]] poisson_interval(
            double b, double e, double m)
            : genetic_map_unit(), beg(b), end(e), mean(m),
              cast(detail::generate_genetic_map_unit_cast_function(false)),
              discrete_(false)
        {
            if (!std::isfinite(beg))
                {
                    throw std::invalid_argument("beg must be finite");
                }
            if (!std::isfinite(end))
                {
                    throw std::invalid_argument("end must be finite");
                }
            if (!std::isfinite(mean))
                {
                    throw std::invalid_argument("mean must be finite");
                }
            if (end <= beg)
                {
                    throw std::invalid_argument("end must be greater than beg");
                }
            if (mean < 0)
                {
                    throw std::invalid_argument("mean must be non-negative");
                }
        }

        poisson_interval(double b, double e, double m, bool discrete)
            : genetic_map_unit(), beg(b), end(e), mean(m),
              cast(detail::generate_genetic_map_unit_cast_function(discrete)),
              discrete_(discrete)
        {
            if (!std::isfinite(beg))
                {
                    throw std::invalid_argument("beg must be finite");
                }
            if (!std::isfinite(end))
                {
                    throw std::invalid_argument("end must be finite");
                }
            if (!std::isfinite(mean))
                {
                    throw std::invalid_argument("mean must be finite");
                }
            if (end <= beg)
                {
                    throw std::invalid_argument("end must be greater than beg");
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
            for (unsigned i = 0; i < n; ++i)
                {
                    breakpoints.push_back(cast(gsl_ran_flat(r, beg, end)));
                }
        }

        std::unique_ptr<genetic_map_unit>
        clone() const final
        {
            return std::make_unique<poisson_interval>(beg, end, mean, discrete_);
        }

        virtual bool
        discrete() const final
        {
            return discrete_;
        }
    };
} // namespace fwdpp

#endif
