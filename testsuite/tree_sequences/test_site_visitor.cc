#include <iostream>
#include <algorithm>
#include <fwdpp/ts/site_visitor.hpp>
#include <boost/test/unit_test.hpp>
#include "simple_table_collection_infinite_sites.hpp"

BOOST_AUTO_TEST_SUITE(test_site_visitor)

BOOST_FIXTURE_TEST_CASE(test_simple_iteration, simple_table_collection_infinite_sites)
{
    std::size_t num_sites = 0;
    fwdpp::ts::site_visitor<fwdpp::ts::std_table_collection> sv(tables, samples);
    decltype(sv()) i;
    while ((i = sv()) != end(sv))
        {
            ++num_sites;
            auto mr = sv.get_mutations();
            BOOST_REQUIRE_EQUAL(std::distance(mr.first, mr.second), 1);
        }
    BOOST_REQUIRE_EQUAL(num_sites, tables.sites.size());
}

BOOST_AUTO_TEST_SUITE_END()

