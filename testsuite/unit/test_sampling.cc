// Test the construction of populations
// from user input data

#include <cassert>
#include <numeric>
#include <boost/test/unit_test.hpp>
#include <fwdpp/diploid_population.hpp>
#include <fwdpp/popgenmut.hpp>
#include <fwdpp/sampling_functions.hpp>
#include <testsuite/fixtures/sugar_fixtures.hpp>
#include <testsuite/util/quick_evolve_sugar.hpp>
#include <gsl/gsl_matrix_char.h>

BOOST_FIXTURE_TEST_SUITE(test_sampling, diploid_population_objects)

BOOST_AUTO_TEST_CASE(test_taking_sample)
{
    diploids.emplace_back(0, 1);
    haploid_genomes[0].n++;
    haploid_genomes[1].n++;
    poptype p(diploids, haploid_genomes, mutations);
    auto m = fwdpp::sample_individuals(p, std::vector<std::size_t>{ 0, 1 },
                                       true, true, true);
    for (auto k : m.neutral_keys)
        {
            BOOST_REQUIRE_EQUAL(p.mutations[k].neutral, true);
        }
    for (auto k : m.selected_keys)
        {
            BOOST_REQUIRE_EQUAL(p.mutations[k].neutral, false);
        }
    BOOST_REQUIRE_EQUAL(m.neutral_keys.size(), 1);
    BOOST_REQUIRE_EQUAL(m.selected_keys.size(), 2);
    BOOST_REQUIRE_EQUAL(m.neutral.data.size(), 4);
    BOOST_REQUIRE_EQUAL(m.selected.data.size(), 8);
    BOOST_REQUIRE_EQUAL(
        std::accumulate(m.neutral.data.begin(), m.neutral.data.end(), 0), 2);
    BOOST_REQUIRE_EQUAL(
        std::accumulate(m.selected.data.begin(), m.selected.data.end(), 0), 4);
    for (std::size_t i = 0; i < m.neutral.data.size(); ++i)
        {
            if (i % 2 == 0.0)
                {
                    BOOST_REQUIRE_EQUAL(m.neutral.data[i], 1);
                }
            else
                {
                    BOOST_REQUIRE_EQUAL(m.neutral.data[i], 0);
                }
        }
}

BOOST_AUTO_TEST_SUITE_END()
