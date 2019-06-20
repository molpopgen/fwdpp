/*
  \file test_simparams.cc
  Tests for fwdpp/simparams.hpp
*/

#include <iostream>
#include <config.h>
#include <boost/test/unit_test.hpp>
#include <fwdpp/simparams.hpp>

//TODO: implement fixtures
#include <fwdpp/diploid_population.hpp>
#include <fwdpp/recbinder.hpp>
#include <fwdpp/poisson_xover.hpp>
#include <fwdpp/fitness_models.hpp>
#include <fwdpp/popgenmut.hpp>
#include <fwdpp/GSLrng_t.hpp>

using poptype = fwdpp::diploid_population<fwdpp::popgenmut>;
using GSLrng = fwdpp::GSLrng_t<fwdpp::GSL_RNG_MT19937>;

BOOST_AUTO_TEST_SUITE(test_simparams)

BOOST_AUTO_TEST_CASE(test_compilation)
{
    poptype pop(1000);
    GSLrng r(42);
    unsigned generation = 0;
    auto mmodel
        = [&pop, &r, &generation](fwdpp::flagged_mutation_queue &recbin,
                                  poptype::mcont_t &mutations) {
              return fwdpp::infsites_popgenmut(
                  recbin, mutations, r.get(), pop.mut_lookup, generation, 0.0,
                  [&r]() { return gsl_rng_uniform(r.get()); },
                  []() { return 0.0; }, []() { return 0.0; });
          };
    const double littler = 10;
    auto rec
        = fwdpp::recbinder(fwdpp::poisson_xover(littler, 0., 1.), r.get());
    auto gv = fwdpp::multiplicative_diploid(fwdpp::fitness(2.));
    auto params = fwdpp::make_genetic_parameters(
        std::move(gv), std::move(mmodel), std::move(rec));
    static_assert(
        fwdpp::traits::is_rec_model<
            decltype(params.generate_breakpoints), typename poptype::diploid_t,
            typename poptype::haploid_genome_t, typename poptype::mcont_t>::value,
        "invalid recombination model");
    static_assert(std::is_same<typename std::remove_const<decltype(
                                   params.interlocus_recombination)>::type,
                               std::nullptr_t>::value,
                  "interlocus recombination type must be std::nullptr_t");
}

BOOST_AUTO_TEST_SUITE_END()
