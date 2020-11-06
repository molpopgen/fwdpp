#ifndef FWDPP_POISSON_POINT_HPP__
#define FWDPP_POISSON_POINT_HPP__

#include <cmath>
#include <stdexcept>
#include <type_traits>
#include <gsl/gsl_randist.h>
#include "genetic_map_unit.hpp"

namespace fwdpp
{
    /// @struct poisson_point_t poisson_point.hpp fwdpp/genetic_map/poisson_point.hpp
    /// @brief Generate a breakpoint at a position if a Poisson deviate is odd.
    /// @ingroup genetic_map_types
    template <typename T> struct poisson_point_t : public genetic_map_unit
    {
        static_assert(std::is_arithmetic<T>::value, "Template type must be numeric");
        /// Breakpoint position
        const T position;
        /// Mean number of crossovers at @a position
        double mean;
        poisson_point_t(const T pos, const double m)
            : genetic_map_unit(), position(pos), mean(m)
        /// \param pos Breakpoint position
        /// \param m The mean of a Poisson distribution
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
                    breakpoints.push_back(static_cast<double>(position));
                }
        }

        std::unique_ptr<genetic_map_unit>
        clone() const final
        /// Clone object
        /// @return ``std::unique_ptr`` to genetic_map_unit.
        {
            return std::unique_ptr<genetic_map_unit>(
                new poisson_point_t(position, mean));
        }
    };

    /// \typedef poisson_point
    /// \brief A poisson_point_t using a double to represent breakpoints.
    /// \ingroup genetic_map_types
    using poisson_point = poisson_point_t<double>;
} // namespace fwdpp

#endif
