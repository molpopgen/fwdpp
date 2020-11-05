#ifndef FWDPP_BINOMIAL_INTERVAL_HPP__
#define FWDPP_BINOMIAL_INTERVAL_HPP__

#include <cmath>
#include <stdexcept>
#include <type_traits>
#include <gsl/gsl_randist.h>
#include "genetic_map_unit.hpp"

namespace fwdpp
{
    template <typename T> struct binomial_interval_t : public genetic_map_unit
    {
        static_assert(std::is_arithmetic<T>::value, "Template type must be numeric");
        double beg, end, prob;
        binomial_interval_t(T b, T e, double p)
            : genetic_map_unit(), beg(static_cast<double>(b)),
              end(static_cast<double>(e)), prob(p)
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
                    throw std::invalid_argument("end must be greater than beg");
                }
            if (prob < 0)
                {
                    throw std::invalid_argument("prob must be non-negative");
                }
        }

        void
        operator()(const gsl_rng* r, std::vector<double>& breakpoints) const final
        {
            if (gsl_rng_uniform(r) <= prob)
                {
                    auto breakpoint = static_cast<T>(gsl_ran_flat(r, beg, end));
                    breakpoints.push_back(static_cast<double>(breakpoint));
                }
        }

        std::unique_ptr<genetic_map_unit>
        clone() const final
        {
            return std::unique_ptr<genetic_map_unit>(
                new binomial_interval_t(static_cast<T>(beg), static_cast<T>(end), prob));
        }
    };

    using binomial_interval = binomial_interval_t<double>;
} // namespace fwdpp

#endif
