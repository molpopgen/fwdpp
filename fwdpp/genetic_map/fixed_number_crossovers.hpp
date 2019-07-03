#ifndef FWDPP_REGIONS_FIXEDCROSSOVERS_HPP
#define FWDPP_REGIONS_FIXEDCROSSOVERS_HPP

#include <cmath>
#include <stdexcept>
#include <gsl/gsl_randist.h>
#include "genetic_map_unit.hpp"

namespace fwdpp
{
    struct fixed_number_crossovers : public genetic_map_unit
    {
        const double beg, end;
        const int nxovers;
        fixed_number_crossovers(double b, double e, int n)
            : beg(b), end(e), nxovers(n)
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
                    throw std::invalid_argument(
                        "end must be greater than beg");
                }
            if (nxovers < 0)
                {
                    throw std::invalid_argument(
                        "number of crossovers must be non-negative");
                }
        }

        void
        operator()(const gsl_rng* r,
                   std::vector<double>& breakpoints) const final
        {
            for (int i = 0; i < nxovers; ++i)
                {
                    breakpoints.push_back(gsl_ran_flat(r, beg, end));
                }
        }

        std::unique_ptr<genetic_map_unit>
        clone() const final
        {
            return std::unique_ptr<genetic_map_unit>(
                new fixed_number_crossovers(beg, end, nxovers));
        }
    };
} // namespace fwdpp

#endif

