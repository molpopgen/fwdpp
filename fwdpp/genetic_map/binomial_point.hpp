#ifndef FWDPP_BINOMIAL_POINT_HPP
#define FWDPP_BINOMIAL_POINT_HPP

#include <cmath>
#include <stdexcept>
#include <type_traits>
#include <gsl/gsl_randist.h>
#include "genetic_map_unit.hpp"

namespace fwdpp
{
    template <typename T> struct binomial_point_t : public genetic_map_unit
    {
        static_assert(std::is_arithmetic<T>::value, "Template type must be numeric");
        T position;
        double prob;
        binomial_point_t(const T pos, const double pr) : position(pos), prob(pr)
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
        {
            return std::unique_ptr<genetic_map_unit>(
                new binomial_point_t(position, prob));
        }
    };

    using binomial_point = binomial_point_t<double>;
} // namespace fwdpp

#endif

