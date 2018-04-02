#ifndef FWDPP_POISSON_XOVER_HPP__
#define FWDPP_POISSON_XOVER_HPP__

#include <vector>
#include <algorithm>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

namespace fwdpp
{
    struct poisson_xover
    /*!
      \brief Simple model of crossing-over.
      Generates a Poisson-distributed number of recombination breakpoints with
      mean \a recrate that
      are uniformly-distributed between \a minpos and \a maxpos
     */
    {
        const double recrate, minpos, maxpos;

        explicit poisson_xover(const double recrate_, const double minpos_,
                               const double maxpos_)
            : recrate{ recrate_ }, minpos{ minpos_ }, maxpos{ maxpos_ }
        /*!
          \param recrate_ The recombination rate (per diploid_, per region)
          \param minpos_ The minimum recombination position allowed
          \param maxpos_ The maximum recombination position allowed
          the gametes & mutations involve in an x-over.
         */
        {
        }

        poisson_xover(const poisson_xover &) = default;
        poisson_xover(poisson_xover &&) = default;

        std::vector<double>
        operator()(const gsl_rng * r) const
        {
            unsigned nbreaks
                = (recrate > 0) ? gsl_ran_poisson(r, recrate) : 0u;
            if (!nbreaks)
                return {};

            std::vector<double> pos;
            pos.reserve(nbreaks + 1);
            for (unsigned i = 0; i < nbreaks; ++i)
                {
                    pos.emplace_back(gsl_ran_flat(r, minpos, maxpos));
                }
            std::sort(pos.begin(), pos.end());
            // Note: this is required for all vectors of breakpoints!
            pos.emplace_back(std::numeric_limits<double>::max());
            return pos;
        }
    };
}
#endif
