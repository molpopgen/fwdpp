#ifndef FWDPP_TESTSUIT_UTIL_MIGPOP_HPP
#define FWDPP_TESTSUIT_UTIL_MIGPOP_HPP

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

inline std::size_t
migpop(const std::size_t &source_pop, const gsl_rng *r, const double &mig_prob)
/*!
  \brief Quick migration model for two-deme simulations
  \ingroup unit
*/
{
    if (gsl_rng_uniform(r) < mig_prob)
        {
            return !source_pop;
        }
    return source_pop;
}

#endif
