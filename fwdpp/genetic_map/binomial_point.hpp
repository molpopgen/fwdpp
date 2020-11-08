#ifndef FWDPP_BINOMIAL_POINT_HPP
#define FWDPP_BINOMIAL_POINT_HPP

#include <cmath>
#include <stdexcept>
#include <type_traits>
#include <gsl/gsl_randist.h>
#include "genetic_map_unit.hpp"

namespace fwdpp
{
    /// @struct binomial_point_t binomial_point.hpp fwdpp/genetic_map//binomial_point.hpp
    /// @brief Generate a breakpoint at a specific position with a given probability
    /// @ingroup genetic_map_types
    template <typename T> struct binomial_point_t : public genetic_map_unit
    {
        static_assert(std::is_arithmetic<T>::value, "Template type must be numeric");
        /// Breakpoint position
        T position;
        /// Breakpoint probability
        double prob;
        binomial_point_t(const T pos, const double pr) : position(pos), prob(pr)
        /// @param pos The breakpoint position
        /// @param pr The probability of a breakpoint
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
                    breakpoints.push_back(static_cast<double>(position));
                }
        }

        std::unique_ptr<genetic_map_unit>
        clone() const final
        /// Clone object
        /// @return ``std::unique_ptr`` to genetic_map_unit.
        {
            return std::unique_ptr<genetic_map_unit>(
                new binomial_point_t(position, prob));
        }
    };

    /// @typedef binomial_point
    /// @brief A fwdpp::binomial_point_t using double to represent breakpoints
    /// @ingroup genetic_map_types
    using binomial_point = binomial_point_t<double>;
} // namespace fwdpp

#endif

