#include <iostream>
#include <boost/test/unit_test.hpp>
#include <fwdpp/ts/std_table_collection.hpp>
#include <fwdpp/ts/table_collection_functions.hpp>

BOOST_AUTO_TEST_SUITE(test_mutation_and_site_tables)

BOOST_AUTO_TEST_CASE(test_invalid_length)
{
    BOOST_REQUIRE_THROW(fwdpp::ts::std_table_collection(0.), std::invalid_argument);
    BOOST_REQUIRE_THROW(fwdpp::ts::std_table_collection(-1.), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(test_rebuild_one_mutation_per_site_inconsistent_data)
{
    fwdpp::ts::std_table_collection tables(1.);
    tables.emplace_back_site(0.1, fwdpp::ts::default_ancestral_state);
    tables.emplace_back_site(0.1, fwdpp::ts::default_ancestral_state);
    // Mutation refers to site 1
    tables.emplace_back_mutation(0, 0lu, 1lu, fwdpp::ts::default_derived_state,
                                 true);
    fwdpp::ts::sort_mutation_table_and_rebuild_site_table(tables);
    BOOST_REQUIRE_EQUAL(tables.sites.size(), 1);
    BOOST_REQUIRE_EQUAL(tables.mutations.size(), 1);
    BOOST_REQUIRE_EQUAL(tables.mutations[0].site, 0);
}

BOOST_AUTO_TEST_CASE(
    test_rebuild_multiple_mutations_per_site_inconsistent_data)
{
    fwdpp::ts::std_table_collection tables(1.);
    auto site
        = tables.emplace_back_site(0.1, fwdpp::ts::default_ancestral_state);
    tables.emplace_back_mutation(0, 0lu, site,
                                 fwdpp::ts::default_derived_state, true);
    site = tables.emplace_back_site(0.1, fwdpp::ts::default_ancestral_state);
    tables.emplace_back_mutation(1, 1lu, site, std::int8_t{ 2 }, true);
    fwdpp::ts::sort_mutation_table_and_rebuild_site_table(tables);
    BOOST_REQUIRE_EQUAL(tables.sites.size(), 1);
    BOOST_REQUIRE_EQUAL(tables.mutations.size(), 2);
    BOOST_REQUIRE_EQUAL(tables.mutations[0].site, 0);
    BOOST_REQUIRE_EQUAL(tables.mutations[1].site, 0);
}

BOOST_AUTO_TEST_SUITE_END()

