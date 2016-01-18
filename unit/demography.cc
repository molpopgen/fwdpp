/*! \file demography.cc
  \ingroup unit 
  \brief Testing functions in fwdpp/demography.hpp
*/
#define BOOST_TEST_MODULE demography
#define BOOST_TEST_DYN_LINK 

#include <config.h>
#include <algorithm>
#include <boost/test/unit_test.hpp>
#include <fwdpp/diploid.hh>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/sugar/metapop.hpp>
#include <fwdpp/sugar/singlepop.hpp>
#include <fwdpp/sugar/infsites.hpp>
#include <fwdpp/sugar/serialization.hpp>

using mtype = KTfwd::popgenmut;
using singlepop_t = KTfwd::singlepop<mtype>;
using metapop_t = KTfwd::metapop<mtype>;


