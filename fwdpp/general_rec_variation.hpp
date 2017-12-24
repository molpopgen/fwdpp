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
    /*!
	 * Define an interval in which
	 * recombination events occur at a
	 * Poisson rate and generate breakpoint
	 * positions uniformly
	 */
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
            if (mean < 0.)
                {
                    throw std::invalid_argument("mean must be non-negative");
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

        template <typename... args>
        inline void
        operator()(std::vector<double>& breakpoints, args&&...) const
        /// Variadic to be combatible with richer recombination policy requirements
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
    /*!
	 * A callable object that returns a recombination
	 * breakpoint at a specific position with a 
	 * given probability
	 */
    {
        const gsl_rng* r;
        const double prob, pos;
        crossover_point(const gsl_rng* r_, const double rate,
                        const double pos_, const double poisson = true)
            /*!
			 * \param r_ A gsl_rng that must be initialized.
			 * \param rate_ The crossover rate. See note below.
			 * \param pos_ The crossover position to return.
			 * \param poisson Whether \a rate represents a Poisson process (true) or a binomial process (false).
			 *
			 * \note If \a poisson is true, it is converted into a probability of 
			 * recombination using Haldane's mapping function.
			 */
            : r{ r_ },
              prob{ (poisson) ? 0.5 * (1. - std::exp(-2.0 * rate)) : rate },
              pos{ pos_ }
        {
            if (rate < 0. || prob < 0.)
                {
                    throw std::invalid_argument(
                        "recombination probability must be non-negative");
                }
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

        template <typename... args>
        inline void
        operator()(std::vector<double>& breakpoints, args&&...) const
        /// Variadic to be combatible with richer recombination policy requirements
        {
            if (gsl_ran_binomial(r, prob, 1))
                {
                    breakpoints.push_back(pos);
                }
        }
    };

    struct general_rec_variation
    /*!
	 * A generalized genetic map.
	 * It holds a vector of functions that add recombination
	 * breakpoints to a vector.  Examples of such functions
	 * are KTfwd::poisson_interval and KTfwd::crossover_point.
	 */
    {
        std::vector<std::function<void(std::vector<double>&)>> recmap;

        general_rec_variation() : recmap{}
        {
        }

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
