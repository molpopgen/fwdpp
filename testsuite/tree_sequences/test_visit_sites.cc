#include <fwdpp/ts/visit_sites.hpp>
#include <boost/test/unit_test.hpp>
#include "simple_table_collection_infinite_sites.hpp"

BOOST_AUTO_TEST_SUITE(test_visit_sites)

BOOST_FIXTURE_TEST_CASE(test_simple_iteration,
                        simple_table_collection_infinite_sites)
{
    std::size_t num_sites = 0;
    auto f =
        [&num_sites](const fwdpp::ts::marginal_tree& /*m*/,
                     const fwdpp::ts::site& /*s*/,
                     const fwdpp::ts::mutation_key_vector::const_iterator b,
                     const fwdpp::ts::mutation_key_vector::const_iterator e) {
            if (std::distance(b, e) != 1)
                {
                    throw std::runtime_error("incorrect number of mutations");
                }
            ++num_sites;
        };
    BOOST_REQUIRE_NO_THROW({
        fwdpp::ts::visit_sites(tables, samples, f, 0, tables.genome_length());
    });
    BOOST_REQUIRE_EQUAL(num_sites, tables.site_table.size());
}

BOOST_AUTO_TEST_SUITE_END()

