// Test the construction of populations
// from user input data

#include <cassert>
#include <boost/test/unit_test.hpp>
#include <fwdpp/diploid_population.hpp>
#include <fwdpp/popgenmut.hpp>
#include <testsuite/fixtures/sugar_fixtures.hpp>

BOOST_AUTO_TEST_SUITE(test_construct_pops)

BOOST_FIXTURE_TEST_CASE(test_diploid_population, diploid_population_objects)
{
    BOOST_REQUIRE_EQUAL(mutations[0].neutral, true);
    BOOST_REQUIRE_EQUAL(mutations[1].neutral, false);
    BOOST_REQUIRE_EQUAL(mutations[2].neutral, false);

    BOOST_REQUIRE_NO_THROW(poptype p(diploids, haploid_genomes, mutations););
}

BOOST_FIXTURE_TEST_CASE(test_diploid_population_N, diploid_population_objects)
{
    poptype p(diploids, haploid_genomes, mutations);
    BOOST_REQUIRE_EQUAL(p.N, p.diploids.size());
}

BOOST_FIXTURE_TEST_CASE(test_diploid_population_bad_haploid_genome_counts,
                        diploid_population_objects)
{
    // Now, make incorrect data to trigger exceptions
    diploids[0].second = 0;
    BOOST_REQUIRE_THROW(poptype p(diploids, haploid_genomes, mutations),
                        std::runtime_error);
}

BOOST_FIXTURE_TEST_CASE(test_diploid_population_bad_haploid_genome_counts2,
                        diploid_population_objects)
{
    // make the haploid_genome data incorrect
    haploid_genomes[0].n++;
    BOOST_REQUIRE_THROW(poptype p(diploids, haploid_genomes, mutations),
                        std::runtime_error);
}

BOOST_FIXTURE_TEST_CASE(test_diploid_population_incorrect_mutation_storage,
                        diploid_population_objects)
{
    std::swap(haploid_genomes[0].mutations, haploid_genomes[0].smutations);
    BOOST_REQUIRE_THROW(poptype p(diploids, haploid_genomes, mutations),
                        std::logic_error);
}

BOOST_AUTO_TEST_SUITE_END()
