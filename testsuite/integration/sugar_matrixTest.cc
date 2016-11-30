/* \brief Integration and unit tests of data matrix generation
 * \ingroup unit
 */

#include <config.h>
#include <boost/test/unit_test.hpp>
#include "../fixtures/sugar_fixtures.hpp"
#include "../util/quick_evolve_sugar.hpp"
#include <fwdpp/sugar/matrix.hpp>

BOOST_AUTO_TEST_CASE(singlepop_hapmatrix)
{
    using spoptype = singlepop_popgenmut_fixture::poptype;
    spoptype pop(1000);
    simulate_singlepop(pop, 1000);
	auto keys = mutation_keys(pop,{0,1,2,3},true,true);
	auto m = haplotype_matrix(pop,{0,1,2,3},keys.first,keys.second);

	BOOST_REQUIRE_EQUAL(m.nrow,8);
	BOOST_REQUIRE_EQUAL(m.neutral.size(),m.nrow*keys.first.size());
	BOOST_REQUIRE_EQUAL(m.selected.size(),m.nrow*keys.second.size());
	BOOST_REQUIRE_EQUAL(keys.first.size(),m.neutral_positions.size());
	BOOST_REQUIRE_EQUAL(keys.first.size(),m.neutral_popfreq.size());
	BOOST_REQUIRE_EQUAL(keys.second.size(),m.selected_positions.size());
	BOOST_REQUIRE_EQUAL(keys.second.size(),m.selected_popfreq.size());
}

BOOST_AUTO_TEST_CASE(singlepop_genotype_matrix)
{
    using spoptype = singlepop_popgenmut_fixture::poptype;
    spoptype pop(1000);
    simulate_singlepop(pop, 1000);
	auto keys = mutation_keys(pop,{0,1,2,3},true,true);
	auto m = genotype_matrix(pop,{0,1,2,3},keys.first,keys.second);

	BOOST_REQUIRE_EQUAL(m.nrow,4);
	BOOST_REQUIRE_EQUAL(m.neutral.size(),m.nrow*keys.first.size());
	BOOST_REQUIRE_EQUAL(m.selected.size(),m.nrow*keys.second.size());
	BOOST_REQUIRE_EQUAL(keys.first.size(),m.neutral_positions.size());
	BOOST_REQUIRE_EQUAL(keys.first.size(),m.neutral_popfreq.size());
	BOOST_REQUIRE_EQUAL(keys.second.size(),m.selected_positions.size());
	BOOST_REQUIRE_EQUAL(keys.second.size(),m.selected_popfreq.size());
}
