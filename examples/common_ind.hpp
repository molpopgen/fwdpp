#ifndef __FWDPP_EXAMPLES_COMMON_IND_HPP__
#define __FWDPP_EXAMPLES_COMMON_IND_HPP__

#include <config.h>
#include <iostream>
#if defined(USE_BOOST_CONTAINERS)
#define FWDPP_SUGAR_USE_BOOST
#endif

/*
   The various examples will define the appropriate symbol,
   so that the minimum stuff is included

   We expect that the typedef "mtype" exists prior
   to including this header, and is an alias to the simulation's
   mutation type
*/
#ifdef DIPLOID_POPULATION_SIM
#include <fwdpp/diploid_population.hpp>
using diploid_population = fwdpp::diploid_population<mtype>;
#endif

// RNG type
#include <fwdpp/GSLrng_t.hpp>

using GSLrng = fwdpp::GSLrng_t<fwdpp::GSL_RNG_MT19937>;

#endif
