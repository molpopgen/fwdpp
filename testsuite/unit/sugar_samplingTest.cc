/*!
  \file sugar_samplingTest.cc
  \ingroup unit
  \brief Testing fwdpp::sample and fwdpp::sample_separate
*/

#include <config.h>
#include <boost/test/unit_test.hpp>
#include <testsuite/fixtures/sugar_fixtures.hpp>
#include <testsuite/util/quick_evolve_sugar.hpp>
#include <fwdpp/sugar/sampling.hpp>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/debug.hpp>
#include <iostream>

fwdpp::GSLrng_t<fwdpp::GSL_RNG_TAUS2> rng(0u);

BOOST_FIXTURE_TEST_SUITE(test_slocuspop_sampling, slocuspop_popgenmut_fixture)

BOOST_AUTO_TEST_CASE(singledeme_test_sep_empty)
{
    auto s = fwdpp::sample_separate(pop, std::vector<unsigned>(), true);
    BOOST_REQUIRE(s.first.empty() == true);
    BOOST_REQUIRE(s.second.empty() == true);
}

BOOST_AUTO_TEST_CASE(singledeme_test_sep_throw)
{
    BOOST_REQUIRE_THROW(auto s = fwdpp::sample_separate(
                            pop, std::vector<unsigned>({ 1000 }), true),
                        std::out_of_range);
}

BOOST_AUTO_TEST_CASE(singledeme_test_empty)
{
    auto s = fwdpp::sample(pop, std::vector<unsigned>(), true);
    BOOST_REQUIRE(s.empty() == true);
}

BOOST_AUTO_TEST_CASE(singledeme_test_throw)
{
    BOOST_REQUIRE_THROW(
        auto s = fwdpp::sample(pop, std::vector<unsigned>({ 1000 }), true),
        std::out_of_range);
}

// Correctness tests

BOOST_AUTO_TEST_CASE(slocuspop_1)
{
    pop.mutations.emplace_back(0.1, 0, 0, 0);
    pop.mcounts.emplace_back(1);
    pop.gametes.emplace_back(
        1,
        std::vector<slocuspop_popgenmut_fixture::poptype::gamete_t::index_t>{
            0 },
        std::
            vector<slocuspop_popgenmut_fixture::poptype::gamete_t::index_t>{});
    pop.gametes[0].n--;
    BOOST_REQUIRE_EQUAL(fwdpp::check_sum(pop.gametes, 2000), true);
    BOOST_REQUIRE_EQUAL(fwdpp::popdata_sane(pop.diploids, pop.gametes,
                                            pop.mutations, pop.mcounts),
                        true);
    pop.diploids[1].first = 1;

    auto x = fwdpp::sample(pop, { 0, 1 }, true);

    BOOST_REQUIRE_EQUAL(x.size(), 1);
    BOOST_REQUIRE_EQUAL(
        std::count(x[0].second.begin(), x[0].second.end(), '1'), 1);
    BOOST_REQUIRE_EQUAL(x[0].second[2], '1');
    // Now, make that diploid homozygous for the mutation
    /*
      FROM HERE ON OUT, INTERNAL DATA STRUCTURES ARE INCONSISTENT.
      This doesn't matter for these tests, but it could in future releases,
      e.g.,
      if internal code adds some assertions...
    */
    pop.diploids[1].second = 1;

    x = fwdpp::sample(pop, { 0, 1 }, true);
    BOOST_REQUIRE_EQUAL(x.size(), 1);
    BOOST_REQUIRE_EQUAL(
        std::count(x[0].second.begin(), x[0].second.end(), '1'), 2);
    BOOST_REQUIRE_EQUAL(x[0].second[2], '1');
    BOOST_REQUIRE_EQUAL(x[0].second[3], '1');

    // Now, we are only sampling individual 1, and thus
    // the variant at 0.1 will be fixed in the sample,
    // and the sample should be empty
    x = fwdpp::sample(pop, { 1 }, true);
    BOOST_REQUIRE_EQUAL(x.empty(), true);

    // Now, allow fixed variants in the sample
    x = fwdpp::sample(pop, { 1 }, false);
    BOOST_REQUIRE_EQUAL(x.empty(), false);

    // now, add a fixation with position < 0.1
    pop.fixations.emplace_back(-0.1, 0, 0, 0);

    // Now, we sample, but we don't want the fixations...
    x = fwdpp::sample(pop, { 0, 1 }, true);
    BOOST_REQUIRE_EQUAL(x.size(), 1);
    BOOST_REQUIRE_EQUAL(
        std::count(x[0].second.begin(), x[0].second.end(), '1'), 2);
    BOOST_REQUIRE_EQUAL(x[0].second[2], '1');
    BOOST_REQUIRE_EQUAL(x[0].second[3], '1');

    // Now, we want the fixation
    x = fwdpp::sample(pop, { 0, 1 }, false);
    BOOST_REQUIRE_EQUAL(x.size(), 2);
    BOOST_REQUIRE_EQUAL(x[0].first, -0.1);
    BOOST_REQUIRE_EQUAL(x[1].first, 0.1);
    BOOST_REQUIRE_EQUAL(std::is_sorted(x.begin(), x.end(),
                                       [](const fwdpp::sample_site_t &i,
                                          const fwdpp::sample_site_t &j) {
                                           return i.first < j.first;
                                       }),
                        true);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(test_slocuspop_custom_sampling,
                         slocuspop_popgenmut_custom_fixture)

BOOST_AUTO_TEST_CASE(singledeme_test_sep_empty)
{
    auto s = fwdpp::sample_separate(pop, std::vector<unsigned>(), true);
    BOOST_REQUIRE(s.first.empty() == true);
    BOOST_REQUIRE(s.second.empty() == true);
}

BOOST_AUTO_TEST_CASE(singledeme_test_sep_throw)
{
    BOOST_REQUIRE_THROW(auto s = fwdpp::sample_separate(
                            pop, std::vector<unsigned>({ 1000 }), true),
                        std::out_of_range);
}

BOOST_AUTO_TEST_CASE(singledeme_test_empty)
{
    auto s = fwdpp::sample(pop, std::vector<unsigned>(), true);
    BOOST_REQUIRE(s.empty() == true);
}

BOOST_AUTO_TEST_CASE(singledeme_test_throw)
{
    BOOST_REQUIRE_THROW(
        auto s = fwdpp::sample(pop, std::vector<unsigned>({ 1000 }), true),
        std::out_of_range);
}

// Correctness tests

BOOST_AUTO_TEST_CASE(slocuspop_1)
{
    pop.mutations.emplace_back(0.1, 0, 0, 0);
    pop.mcounts.emplace_back(1);
    pop.gametes.emplace_back(
        1,
        std::vector<slocuspop_popgenmut_fixture::poptype::gamete_t::index_t>{
            0 },
        std::
            vector<slocuspop_popgenmut_fixture::poptype::gamete_t::index_t>{});
    pop.gametes[0].n--;
    BOOST_REQUIRE_EQUAL(fwdpp::check_sum(pop.gametes, 2000), true);
    BOOST_REQUIRE_EQUAL(fwdpp::popdata_sane(pop.diploids, pop.gametes,
                                            pop.mutations, pop.mcounts),
                        true);
    pop.diploids[1].first = 1;

    auto x = fwdpp::sample(pop, { 0, 1 }, true);

    BOOST_REQUIRE_EQUAL(x.size(), 1);
    BOOST_REQUIRE_EQUAL(
        std::count(x[0].second.begin(), x[0].second.end(), '1'), 1);
    BOOST_REQUIRE_EQUAL(x[0].second[2], '1');
    // Now, make that diploid homozygous for the mutation
    /*
      FROM HERE ON OUT, INTERNAL DATA STRUCTURES ARE INCONSISTENT.
      This doesn't matter for these tests, but it could in future releases,
      e.g.,
      if internal code adds some assertions...
    */
    pop.diploids[1].second = 1;

    x = fwdpp::sample(pop, { 0, 1 }, true);
    BOOST_REQUIRE_EQUAL(x.size(), 1);
    BOOST_REQUIRE_EQUAL(
        std::count(x[0].second.begin(), x[0].second.end(), '1'), 2);
    BOOST_REQUIRE_EQUAL(x[0].second[2], '1');
    BOOST_REQUIRE_EQUAL(x[0].second[3], '1');

    // Now, we are only sampling individual 1, and thus
    // the variant at 0.1 will be fixed in the sample,
    // and the sample should be empty
    x = fwdpp::sample(pop, { 1 }, true);
    BOOST_REQUIRE_EQUAL(x.empty(), true);

    // Now, allow fixed variants in the sample
    x = fwdpp::sample(pop, { 1 }, false);
    BOOST_REQUIRE_EQUAL(x.empty(), false);

    // now, add a fixation with position < 0.1
    pop.fixations.emplace_back(-0.1, 0, 0, 0);

    // Now, we sample, but we don't want the fixations...
    x = fwdpp::sample(pop, { 0, 1 }, true);
    BOOST_REQUIRE_EQUAL(x.size(), 1);
    BOOST_REQUIRE_EQUAL(
        std::count(x[0].second.begin(), x[0].second.end(), '1'), 2);
    BOOST_REQUIRE_EQUAL(x[0].second[2], '1');
    BOOST_REQUIRE_EQUAL(x[0].second[3], '1');

    // Now, we want the fixation
    x = fwdpp::sample(pop, { 0, 1 }, false);
    BOOST_REQUIRE_EQUAL(x.size(), 2);
    BOOST_REQUIRE_EQUAL(x[0].first, -0.1);
    BOOST_REQUIRE_EQUAL(x[1].first, 0.1);
    BOOST_REQUIRE_EQUAL(std::is_sorted(x.begin(), x.end(),
                                       [](const fwdpp::sample_site_t &i,
                                          const fwdpp::sample_site_t &j) {
                                           return i.first < j.first;
                                       }),
                        true);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(test_multilocus_sampling, mlocuspop_popgenmut_fixture)

BOOST_AUTO_TEST_CASE(multilocus_test_sep_empty)
{
    auto s = fwdpp::sample_separate(rng.get(), pop, 20, true);
    BOOST_REQUIRE_EQUAL(s.size(), 4);
    for (auto i : s)
        {
            BOOST_REQUIRE(i.first.empty() == true);
            BOOST_REQUIRE(i.second.empty() == true);
        }
}

BOOST_AUTO_TEST_CASE(multilocus_test_sampling)
{
	std::cout << recmodels.size() << '\n';
	for(std::size_t i=0;i<recmodels.size();++i)
	{
		std::cout << i << std::endl;
		auto x = recmodels.at(2)();
		std::cout << x.size() << std::endl;
	}
    //simulate_mlocuspop(pop, rng, mutmodels, recmodels, multilocus_additive(),
    //                   mu, rbw, generation);
    //auto s = fwdpp::sample_separate(rng.get(), pop, 20, true);
    //s = fwdpp::sample_separate(rng.get(), pop, 20, false);
    //auto s2 = fwdpp::sample(rng.get(), pop, 20, true);
    //s2 = fwdpp::sample(rng.get(), pop, 20, false);
    //pop.locus_boundaries.clear();
    //BOOST_REQUIRE_THROW(s = fwdpp::sample_separate(rng.get(), pop, 20, false),
    //                    std::runtime_error);
    //BOOST_REQUIRE_THROW(s2 = fwdpp::sample(rng.get(), pop, 20, false),
    //                    std::runtime_error);
}

BOOST_AUTO_TEST_CASE(multilocus_test_empty)
{
    auto s = fwdpp::sample(rng.get(), pop, 20, true);
    BOOST_REQUIRE_EQUAL(s.size(), 4);
    for (auto i : s)
        {
            BOOST_REQUIRE(i.empty() == true);
        }
}

BOOST_AUTO_TEST_SUITE_END()

