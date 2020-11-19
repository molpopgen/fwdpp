/*
  \file extensions_regionsTest.cc
  API checks on fwdpp's extensions sub-library.
*/

#include <config.h>
#include <boost/test/unit_test.hpp>
#include <fwdpp/GSLrng_t.hpp>
#include <fwdpp/extensions/regions.hpp>
#include <fwdpp/type_traits.hpp>
#include <limits>
#include "../fixtures/sugar_fixtures.hpp"

using namespace fwdpp;

BOOST_FIXTURE_TEST_SUITE(test_extensions, diploid_population_mutation_fixture)

// The next two test cases use a simple dependency injection to test concepts.
// The function types defined have the correct signature (args + return value),
// but do not do the correct operations within, which is fine for the purposes
// here.

BOOST_AUTO_TEST_CASE(bind_mutation_model)
{
    using poptype = diploid_population_mutation_fixture::poptype;
    using mmodel_type = fwdpp::traits::mutation_model<poptype::mutation_container>;
    rng_t rng{42};
    const auto mmodel = [](fwdpp::flagged_mutation_queue &,
                           poptype::mutation_container &) -> std::size_t { return 1; };

    fwdpp::extensions::discrete_mut_model<poptype::mutation_container> m(
        std::vector<mmodel_type>{mmodel}, std::vector<double>{1});

    auto bound_mmodel = fwdpp::extensions::bind_dmm(rng.get(), m);

    auto is_a_mutation_model
        = fwdpp::traits::is_mutation_model<decltype(bound_mmodel),
                                           poptype::mutation_container,
                                           poptype::genome_container>::value;

    BOOST_REQUIRE(is_a_mutation_model);
}

BOOST_AUTO_TEST_CASE(bind_haploid_genome_dependent_mutation_model)
{
    using poptype = diploid_population_mutation_fixture::poptype;
    using mmodel_type
        = fwdpp::traits::mutation_model_haploid_genome<poptype::mutation_container,
                                                       poptype::genome_container>;
    rng_t rng{42};
    const auto mmodel
        = [](const poptype::haploid_genome_type &, fwdpp::flagged_mutation_queue &,
             poptype::mutation_container &) -> std::size_t { return 1; };

    fwdpp::extensions::discrete_mut_model<poptype::mutation_container,
                                          poptype::genome_container>
        m(std::vector<mmodel_type>{mmodel}, std::vector<double>{1});

    auto bound_mmodel = fwdpp::extensions::bind_dmm(rng.get(), m);

    auto is_a_mutation_model
        = fwdpp::traits::is_mutation_model<decltype(bound_mmodel),
                                           poptype::mutation_container,
                                           poptype::genome_container>::value;

    BOOST_REQUIRE(is_a_mutation_model);
}

BOOST_AUTO_TEST_SUITE_END()
