#ifndef FWDPP_GENERAL_REC_VARIATION_HPP__
#define FWDPP_GENERAL_REC_VARIATION_HPP__

#include <vector>
#include <functional>
#include <algorithm>
#include <limits>
#include <cmath>
#include <stdexcept>
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
            if (r == nullptr)
                {
                    throw std::invalid_argument("rng cannot be null pointer");
                }
            if (!std::isfinite(mean))
                {
                    throw std::invalid_argument("mean must be finite");
                }
            if (!std::isfinite(minpos))
                {
                    throw std::invalid_argument("minpos must be finite");
                }
            if (!std::isfinite(maxpos))
                {
                    throw std::invalid_argument("maxpos must be finite");
                }
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

    struct crossover_point
    {
        const gsl_rng* r;
        const double prob, pos;
        crossover_point(const gsl_rng* r_, const double rate,
                        const double pos_, const double poisson = true)
            : r{ r_ },
              prob{ (poisson) ? 0.5 * (1. - std::exp(-2.0 * rate)) : rate },
              pos{ pos_ }
        {
            if (r == nullptr)
                {
                    throw std::invalid_argument("rng cannot be null pointer");
                }
            if (!std::isfinite(rate))
                {
                    throw std::invalid_argument("rate must be finite");
                }
            if (!std::isfinite(prob))
                {
                    throw std::invalid_argument("prob must be finite");
                }
            if (!std::isfinite(pos))
                {
                    throw std::invalid_argument("pos must be finite");
                }
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
            std::sort(breakpoints.begin(), breakpoints.end());
            breakpoints.push_back(std::numeric_limits<double>::max());
            return breakpoints;
        }
    };
}

#endif
