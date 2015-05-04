#ifndef __FWDPP_EXAMPLES_COMMON_IND_HPP__
#define __FWDPP_EXAMPLES_COMMON_IND_HPP__

#include <config.h>
#include <iostream>

#if defined(HAVE_BOOST_VECTOR) && defined(HAVE_BOOST_LIST) && defined(HAVE_BOOST_UNORDERED_SET) && defined(HAVE_BOOST_POOL_ALLOC) && defined(HAVE_BOOST_HASH) && !defined(USE_STANDARD_CONTAINERS)
#define FWDPP_SUGAR_USE_BOOST
#endif

/* 
   The various examples will define the appropriate symbol,
   so that the minimum stuff is included

   We expect that the typedef "mtype" exists prior
   to including this header, and is an alias to the simulation's
   mutation type
*/
#include <fwdpp/sugar/serialization.hpp>
#ifdef SINGLEPOP_SIM
#include <fwdpp/sugar/singlepop.hpp>
using singlepop_t = KTfwd::singlepop<mtype>;
using singlepop_serialized_t = KTfwd::singlepop_serialized<mtype,KTfwd::mutation_writer,KTfwd::mutation_reader<mtype>>;
#elif defined(METAPOP_SIM)
#include <fwdpp/sugar/metapop.hpp>
using metapop_t = KTfwd::metapop<mtype>;
using metapop_serialized_t = KTfwd::metapop_serialized<mtype,KTfwd::mutation_writer,KTfwd::mutation_reader<mtype>>;
#elif defined(MULTILOCUS_SIM)
#include <fwdpp/sugar/multiloc.hpp>
using multiloc_t = KTfwd::multiloc<mtype>;
using multiloc_serialized_t = KTfwd::multiloc_serialized<mtype,KTfwd::mutation_writer,KTfwd::mutation_reader<mtype>>;
#endif

//RNG type
#include <fwdpp/sugar/GSLrng_t.hpp>

using GSLrng = KTfwd::GSLrng_t<KTfwd::GSL_RNG_MT19937>;

#endif
