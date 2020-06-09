#include <iostream>
#include <limits>
#include <boost/test/unit_test.hpp>
#include <fwdpp/ts/recording/diploid_offspring.hpp>
#include <fwdpp/ts/std_table_collection.hpp>

BOOST_AUTO_TEST_SUITE(test_diploid_recording)

BOOST_AUTO_TEST_CASE(test_first_breakpoint_at_zero)
{
	std::vector<double> breakpoints{0, 0.5, std::numeric_limits<double>::max()};
	fwdpp::ts::std_table_collection tables{1.};
	tables.emplace_back_node(0, 0.);

	fwdpp::ts::record_diploid_offspring(breakpoints, std::make_tuple(0, 1),
			0, 1., tables);
	BOOST_REQUIRE_EQUAL(tables.num_edges(), 2);
	BOOST_REQUIRE_EQUAL(tables.edges[0].parent, 1);
	BOOST_REQUIRE_EQUAL(tables.edges[0].left, 0);
	BOOST_REQUIRE_EQUAL(tables.edges[1].parent, 0);
	BOOST_REQUIRE_EQUAL(tables.edges[1].left, 0.5);
}

BOOST_AUTO_TEST_SUITE_END()

