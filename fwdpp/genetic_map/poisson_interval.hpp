#ifndef FWDPP_POISSON_INTERVAL_HPP__
#define FWDPP_POISSON_INTERVAL_HPP__

#include <cmath>
#include <stdexcept>
#include <type_traits>
#include <gsl/gsl_randist.h>
#include "genetic_map_unit.hpp"

namespace fwdpp
{
    /// @struct poisson_interval_t poisson_interval.hpp fwdpp/genetic_map/poisson_interval.hpp
    /// @brief Generate a Poisson number of breakpoints in an interval
    /// @ingroup genetic_map_types
    template <typename T> struct poisson_interval_t : public genetic_map_unit
    {
        static_assert(std::is_arithmetic<T>::value, "Template type must be numeric");
        /// Beginning of interval (inclusive)
        const double beg;
        /// End of interval (exclusive)
        const double end;
        /// Mean number of breakpoints
        const double mean;
        poisson_interval_t(T b, T e, double m)
            : genetic_map_unit(), beg(static_cast<double>(b)),
              end(static_cast<double>(e)), mean(m)
        /// @param b Beginning of interval
        /// @param e End of interval
        /// @param m Mean number of breakpoints
        /// @note The interval is half-open on [b, e).
        {
            internal::validate_minimum_interval_size(b, e);
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
                    auto breakpoint = static_cast<T>(gsl_ran_flat(r, beg, end));
                    breakpoints.push_back(static_cast<double>(breakpoint));
                }
        }

        std::unique_ptr<genetic_map_unit>
        clone() const final
        /// Clone object
        /// @return ``std::unique_ptr`` to genetic_map_unit.
        {
            return std::unique_ptr<genetic_map_unit>(
                new poisson_interval_t(static_cast<T>(beg), static_cast<T>(end), mean));
        }
    };

    /// @typedef poisson_interval
    /// @brief A fwdpp::poisson_interval_t using double to represent breakpoints
    /// @ingroup genetic_map_types
    using poisson_interval = poisson_interval_t<double>;
} // namespace fwdpp

#endif
