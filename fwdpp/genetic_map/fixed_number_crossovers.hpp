#ifndef FWDPP_REGIONS_FIXEDCROSSOVERS_HPP
#define FWDPP_REGIONS_FIXEDCROSSOVERS_HPP

#include <cmath>
#include <stdexcept>
#include <type_traits>
#include <gsl/gsl_randist.h>
#include "genetic_map_unit.hpp"

namespace fwdpp
{
    template <typename T> struct fixed_number_crossovers_t : public genetic_map_unit
    {
        static_assert(std::is_arithmetic<T>::value, "Template type must be numeric");
        const double beg, end;
        const int nxovers;
        fixed_number_crossovers_t(T b, T e, int n)
            : beg(static_cast<double>(b)), end(static_cast<double>(e)), nxovers(n)
        {
            if (!std::isfinite(b))
                {
                    throw std::invalid_argument("beg must be finite");
                }
            if (!std::isfinite(e))
                {
                    throw std::invalid_argument("end must be finite");
                }
            if (e <= b)
                {
                    throw std::invalid_argument("end must be greater than beg");
                }
            if (nxovers < 0)
                {
                    throw std::invalid_argument(
                        "number of crossovers must be non-negative");
                }
        }

        void
        operator()(const gsl_rng* r, std::vector<double>& breakpoints) const final
        {
            for (int i = 0; i < nxovers; ++i)
                {
                    auto breakpoint = static_cast<T>(gsl_ran_flat(r, beg, end));
                    breakpoints.push_back(static_cast<double>(breakpoint));
                }
        }

        std::unique_ptr<genetic_map_unit>
        clone() const final
        {
            return std::unique_ptr<genetic_map_unit>(new fixed_number_crossovers_t(
                static_cast<T>(beg), static_cast<T>(end), nxovers));
        }
    };

    using fixed_number_crossovers = fixed_number_crossovers_t<double>;
} // namespace fwdpp

#endif

