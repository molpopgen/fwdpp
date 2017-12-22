#ifndef FWDPP_GENERAL_REC_VARIATION_HPP__
#define FWDPP_GENERAL_REC_VARIATION_HPP__

#include <vector>
#include <functional>
#include <algorithm>
#include <limits>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

namespace KTfwd
{
    struct poisson_interval
    {
        const gsl_rng* r;
        const double mean, minpos, maxpos;
        poisson_interval(const gsl_rng* r_, const double mean_,
                         const double minpos_, const double maxpos_)
            : r{ r_ }, mean{ mean_ }, minpos{ minpos_ }, maxpos{ maxpos_ }
        {
        }

        inline void
        operator()(std::vector<double>& breakpoints) const
        {
            auto nbreaks = gsl_ran_poisson(r, mean);
            for (decltype(nbreaks) i = 0; i < nbreaks; ++i)
                {
                    double pos = gsl_ran_flat(r, minpos, maxpos);
                    breakpoints.push_back(pos);
                }
        }
    };

    struct poisson_point
    {
        const double prob, pos;
        poisson_point(const double rate, const double pos_)
            : prob{ 0.5 * (1. - std::exp(-2.0 * rate)) }, pos{ pos_ }
        {
        }

        inline void
        operator()(std::vector<double>& breakpoints) const
        {
            if (gsl_ran_binomial(r, prob, 1))
                {
                    breakpoints.push_back(pos);
                }
        }
    };

    struct binomial_point
    {
        const double prob, pos;
        poisson_point(const double prob_, const double pos_)
            : prob{ prob_ }, pos{ pos_ }
        {
        }

        inline void
        operator()(std::vector<double>& breakpoints) const
        {
            if (gsl_ran_binomial(r, prob, 1))
                {
                    breakpoints.push_back(pos);
                }
        }
    };

    struct general_rec_variation
    {
        std::vector<std::function<void(std::vector<double>&)>> recmap;

        inline std::vector<double>
        operator()() const
        {
            std::vector<double> breakpoints;
            for (const auto& f : recmap)
                {
                    f(breakpoints);
                }
            std::sort(f.begin(), f.end());
            breakpoints.push_back(std::numeric_limits<double>::max());
            return breakpoints;
        }
    };
}

#endif
