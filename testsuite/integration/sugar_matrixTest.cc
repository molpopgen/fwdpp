/* \brief Integration and unit tests of data matrix generation
 * \ingroup unit
 */

#include <config.h>
#include <boost/test/unit_test.hpp>
#include "../fixtures/sugar_fixtures.hpp"
#include "../util/quick_evolve_sugar.hpp"
#include <fwdpp/sugar/matrix.hpp>
#include <fwdpp/sugar/sampling.hpp>

BOOST_AUTO_TEST_CASE(singlepop_hapmatrix)
{
    using spoptype = singlepop_popgenmut_fixture::poptype;
    spoptype pop(1000);
    simulate_singlepop(pop, 10000);
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

BOOST_AUTO_TEST_CASE(singlepop_hapmatrix_compare_to_sample)
{
    using spoptype = singlepop_popgenmut_fixture::poptype;
    spoptype pop(1000);
    simulate_singlepop(pop, 10000);
	std::vector<std::size_t> indlist;
	//Sample a LOT of individuals
	for( std::size_t i = 100 ; i < 750 ; i += 5 ) indlist.push_back(i);
	std::vector<unsigned> indlist2(indlist.begin(),indlist.end());
	auto keys = mutation_keys(pop,indlist,true,true);
	auto m = haplotype_matrix(pop,indlist,keys.first,keys.second);
	auto s = KTfwd::sample_separate(pop,indlist2,true);
	BOOST_REQUIRE_EQUAL(s.first.size(),m.neutral_positions.size());
	BOOST_REQUIRE_EQUAL(s.second.size(),m.selected_positions.size());
	std::vector<double> pos;
	for(auto && i : s.first) pos.push_back(i.first);
	for(std::size_t i = 0 ; i < pos.size() ; ++i)
	{
		BOOST_REQUIRE_EQUAL(pos[i],m.neutral_positions[i]);
	}
	pos.clear();
	for(auto && i : s.second) pos.push_back(i.first);
	for(std::size_t i = 0 ; i < pos.size() ; ++i)
	{
		BOOST_REQUIRE_EQUAL(pos[i],m.selected_positions[i]);
	}
	BOOST_REQUIRE_EQUAL(m.nrow,2*indlist.size());
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
    simulate_singlepop(pop, 10000);
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
