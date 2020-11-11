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
#include <fwdpp/genetic_map/genetic_map.hpp>
#include <fwdpp/genetic_map/poisson_interval.hpp>
#include <fwdpp/fitness_models.hpp>
#include <fwdpp/types/mutation.hpp>
#include <fwdpp/GSLrng_t.hpp>

using poptype = fwdpp::diploid_population<fwdpp::mutation>;
using GSLrng = fwdpp::GSLrng_t<fwdpp::GSL_RNG_MT19937>;

BOOST_AUTO_TEST_SUITE(test_simparams)

BOOST_AUTO_TEST_CASE(test_compilation)
{
    poptype pop(1000);
    GSLrng r(42);
    unsigned generation = 0;
    auto mmodel
        = [&pop, &r, &generation](fwdpp::flagged_mutation_queue &recbin,
                                  poptype::mutation_container &mutations) {
              return fwdpp::infsites_mutation(
                  recbin, mutations, r.get(), pop.mut_lookup, generation, 0.0,
                  [&r]() { return gsl_rng_uniform(r.get()); },
                  []() { return 0.0; }, []() { return 0.0; });
          };
    const double littler = 10;
    fwdpp::genetic_map gmap;
    gmap.add_callback(fwdpp::poisson_interval(0, 1, littler));
    const auto rec = fwdpp::recbinder(std::cref(gmap), r.get());

    auto gv = fwdpp::multiplicative_diploid(fwdpp::fitness(2.));
    auto params = fwdpp::make_genetic_parameters(
        std::move(gv), std::move(mmodel), std::move(rec));
    static_assert(
        fwdpp::traits::is_rec_model<
            decltype(params.generate_breakpoints), typename poptype::diploid_type,
            typename poptype::haploid_genome_type, typename poptype::mutation_container>::value,
        "invalid recombination model");
}

BOOST_AUTO_TEST_SUITE_END()
