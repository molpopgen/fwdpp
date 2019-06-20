/*!
  \file test_recombination.cc
  \ingroup unit
*/
#include <config.h>
#include <iostream>
// For this unit test, this symbol eliminates the mutation-related part of
// fwdpp::fwdpp_internal::multiloc_rec_mut,
// which means we don't have to write as much boilerplate code to test the more
// complex logic.
// Plus, mutation stuff is unit-tested elsewhere
// #define FWDPP_UNIT_TESTING
#include <fwdpp/diploid.hh>
#include <boost/test/unit_test.hpp>
#include <unistd.h>
#include <iterator>
#include <functional>
#include <vector>
#include <gsl/gsl_rng.h>
#include <fwdpp/mutate_recombine.hpp>

BOOST_AUTO_TEST_SUITE(test_recombination)

struct simple_mutation : public fwdpp::mutation_base
{
    simple_mutation(double p) : fwdpp::mutation_base(p) {}
};


BOOST_AUTO_TEST_SUITE_END()
