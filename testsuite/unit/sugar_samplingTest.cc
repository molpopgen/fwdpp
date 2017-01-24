/*!
  \file sugar_sampling.cc
  \ingroup unit
  \brief Testing KTfwd::sample and KTfwd::sample_separate
*/
#include <config.h>
#include <boost/test/unit_test.hpp>
#include <testsuite/fixtures/sugar_fixtures.hpp>
#include <fwdpp/sugar/sampling.hpp>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/debug.hpp>

KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);

BOOST_FIXTURE_TEST_SUITE(test_singlepop_sampling, singlepop_popgenmut_fixture)

BOOST_AUTO_TEST_CASE(singledeme_test_sep_empty)
{
    auto s = KTfwd::sample_separate(pop, std::vector<unsigned>(), true);
    BOOST_REQUIRE(s.first.empty() == true);
    BOOST_REQUIRE(s.second.empty() == true);
}

BOOST_AUTO_TEST_CASE(singledeme_test_sep_throw)
{
    BOOST_REQUIRE_THROW(auto s = KTfwd::sample_separate(
                            pop, std::vector<unsigned>({ 1000 }), true),
                        std::out_of_range);
}

BOOST_AUTO_TEST_CASE(singledeme_test_empty)
{
    auto s = KTfwd::sample(pop, std::vector<unsigned>(), true);
    BOOST_REQUIRE(s.empty() == true);
}

BOOST_AUTO_TEST_CASE(singledeme_test_throw)
{
    BOOST_REQUIRE_THROW(
        auto s = KTfwd::sample(pop, std::vector<unsigned>({ 1000 }), true),
        std::out_of_range);
}

// Correctness tests

BOOST_AUTO_TEST_CASE(singlepop_1)
{
    pop.mutations.emplace_back(0.1, 0, 0, 0);
    pop.mcounts.emplace_back(1);
    pop.gametes.emplace_back(
        1,
        std::vector<singlepop_popgenmut_fixture::poptype::gamete_t::index_t>{
            0 },
        std::
            vector<singlepop_popgenmut_fixture::poptype::gamete_t::index_t>{});
    pop.gametes[0].n--;
    BOOST_REQUIRE_EQUAL(KTfwd::check_sum(pop.gametes, 2000), true);
    BOOST_REQUIRE_EQUAL(KTfwd::popdata_sane(pop.diploids, pop.gametes,
                                            pop.mutations, pop.mcounts),
                        true);
    pop.diploids[1].first = 1;

    auto x = KTfwd::sample(pop, { 0, 1 }, true);

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

    x = KTfwd::sample(pop, { 0, 1 }, true);
    BOOST_REQUIRE_EQUAL(x.size(), 1);
    BOOST_REQUIRE_EQUAL(
        std::count(x[0].second.begin(), x[0].second.end(), '1'), 2);
    BOOST_REQUIRE_EQUAL(x[0].second[2], '1');
    BOOST_REQUIRE_EQUAL(x[0].second[3], '1');

    // Now, we are only sampling individual 1, and thus
    // the variant at 0.1 will be fixed in the sample,
    // and the sample should be empty
    x = KTfwd::sample(pop, { 1 }, true);
    BOOST_REQUIRE_EQUAL(x.empty(), true);

    // Now, allow fixed variants in the sample
    x = KTfwd::sample(pop, { 1 }, false);
    BOOST_REQUIRE_EQUAL(x.empty(), false);

    // now, add a fixation with position < 0.1
    pop.fixations.emplace_back(-0.1, 0, 0, 0);

    // Now, we sample, but we don't want the fixations...
    x = KTfwd::sample(pop, { 0, 1 }, true);
    BOOST_REQUIRE_EQUAL(x.size(), 1);
    BOOST_REQUIRE_EQUAL(
        std::count(x[0].second.begin(), x[0].second.end(), '1'), 2);
    BOOST_REQUIRE_EQUAL(x[0].second[2], '1');
    BOOST_REQUIRE_EQUAL(x[0].second[3], '1');

    // Now, we want the fixation
    x = KTfwd::sample(pop, { 0, 1 }, false);
    BOOST_REQUIRE_EQUAL(x.size(), 2);
    BOOST_REQUIRE_EQUAL(x[0].first, -0.1);
    BOOST_REQUIRE_EQUAL(x[1].first, 0.1);
    BOOST_REQUIRE_EQUAL(std::is_sorted(x.begin(), x.end(),
                                       [](const KTfwd::sample_site_t &i,
                                          const KTfwd::sample_site_t &j) {
                                           return i.first < j.first;
                                       }),
                        true);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(test_singlepop_custom_sampling,
                         singlepop_popgenmut_custom_fixture)

BOOST_AUTO_TEST_CASE(singledeme_test_sep_empty)
{
    auto s = KTfwd::sample_separate(pop, std::vector<unsigned>(), true);
    BOOST_REQUIRE(s.first.empty() == true);
    BOOST_REQUIRE(s.second.empty() == true);
}

BOOST_AUTO_TEST_CASE(singledeme_test_sep_throw)
{
    BOOST_REQUIRE_THROW(auto s = KTfwd::sample_separate(
                            pop, std::vector<unsigned>({ 1000 }), true),
                        std::out_of_range);
}

BOOST_AUTO_TEST_CASE(singledeme_test_empty)
{
    auto s = KTfwd::sample(pop, std::vector<unsigned>(), true);
    BOOST_REQUIRE(s.empty() == true);
}

BOOST_AUTO_TEST_CASE(singledeme_test_throw)
{
    BOOST_REQUIRE_THROW(
        auto s = KTfwd::sample(pop, std::vector<unsigned>({ 1000 }), true),
        std::out_of_range);
}

// Correctness tests

BOOST_AUTO_TEST_CASE(singlepop_1)
{
    pop.mutations.emplace_back(0.1, 0, 0, 0);
    pop.mcounts.emplace_back(1);
    pop.gametes.emplace_back(
        1,
        std::vector<singlepop_popgenmut_fixture::poptype::gamete_t::index_t>{
            0 },
        std::
            vector<singlepop_popgenmut_fixture::poptype::gamete_t::index_t>{});
    pop.gametes[0].n--;
    BOOST_REQUIRE_EQUAL(KTfwd::check_sum(pop.gametes, 2000), true);
    BOOST_REQUIRE_EQUAL(KTfwd::popdata_sane(pop.diploids, pop.gametes,
                                            pop.mutations, pop.mcounts),
                        true);
    pop.diploids[1].first = 1;

    auto x = KTfwd::sample(pop, { 0, 1 }, true);

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

    x = KTfwd::sample(pop, { 0, 1 }, true);
    BOOST_REQUIRE_EQUAL(x.size(), 1);
    BOOST_REQUIRE_EQUAL(
        std::count(x[0].second.begin(), x[0].second.end(), '1'), 2);
    BOOST_REQUIRE_EQUAL(x[0].second[2], '1');
    BOOST_REQUIRE_EQUAL(x[0].second[3], '1');

    // Now, we are only sampling individual 1, and thus
    // the variant at 0.1 will be fixed in the sample,
    // and the sample should be empty
    x = KTfwd::sample(pop, { 1 }, true);
    BOOST_REQUIRE_EQUAL(x.empty(), true);

    // Now, allow fixed variants in the sample
    x = KTfwd::sample(pop, { 1 }, false);
    BOOST_REQUIRE_EQUAL(x.empty(), false);

    // now, add a fixation with position < 0.1
    pop.fixations.emplace_back(-0.1, 0, 0, 0);

    // Now, we sample, but we don't want the fixations...
    x = KTfwd::sample(pop, { 0, 1 }, true);
    BOOST_REQUIRE_EQUAL(x.size(), 1);
    BOOST_REQUIRE_EQUAL(
        std::count(x[0].second.begin(), x[0].second.end(), '1'), 2);
    BOOST_REQUIRE_EQUAL(x[0].second[2], '1');
    BOOST_REQUIRE_EQUAL(x[0].second[3], '1');

    // Now, we want the fixation
    x = KTfwd::sample(pop, { 0, 1 }, false);
    BOOST_REQUIRE_EQUAL(x.size(), 2);
    BOOST_REQUIRE_EQUAL(x[0].first, -0.1);
    BOOST_REQUIRE_EQUAL(x[1].first, 0.1);
    BOOST_REQUIRE_EQUAL(std::is_sorted(x.begin(), x.end(),
                                       [](const KTfwd::sample_site_t &i,
                                          const KTfwd::sample_site_t &j) {
                                           return i.first < j.first;
                                       }),
                        true);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(test_multilocus_sampling, multiloc_popgenmut_fixture)

BOOST_AUTO_TEST_CASE(multilocus_test_sep_empty)
{
    auto s = KTfwd::sample_separate(rng.get(), pop, 20, true);
    BOOST_REQUIRE_EQUAL(s.size(), 4);
    for (auto i : s)
        {
            BOOST_REQUIRE(i.first.empty() == true);
            BOOST_REQUIRE(i.second.empty() == true);
        }
}

BOOST_AUTO_TEST_CASE(multilocus_test_empty)
{
    auto s = KTfwd::sample(rng.get(), pop, 20, true);
    BOOST_REQUIRE_EQUAL(s.size(), 4);
    for (auto i : s)
        {
            BOOST_REQUIRE(i.empty() == true);
        }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(test_metapop_sampling, metapop_popgenmut_fixture)

BOOST_AUTO_TEST_CASE(metapop_test_sep_deme_out_of_range_1)
{
    BOOST_REQUIRE_THROW(auto s
                        = KTfwd::sample_separate(rng.get(), pop, 2, 10, true),
                        std::out_of_range);
}

BOOST_AUTO_TEST_CASE(metapop_test_sep_deme_out_of_range_2)
{
    BOOST_REQUIRE_THROW(
        auto s = KTfwd::sample_separate(
            pop, 2, std::vector<unsigned>({ 10, 20, 50 }), true),
        std::out_of_range);
}

BOOST_AUTO_TEST_CASE(metapop_test_sep_throw)
{
    BOOST_REQUIRE_THROW(
        auto s = KTfwd::sample_separate(
            pop, 0, std::vector<unsigned>({ 10, 20, 50, 1000 }), true),
        std::out_of_range);
}

BOOST_AUTO_TEST_CASE(metapop_test_sep_empty)
{
    auto s = KTfwd::sample_separate(pop, 0, std::vector<unsigned>(), true);
    BOOST_REQUIRE(s.first.empty() == true);
    BOOST_REQUIRE(s.second.empty() == true);
}

BOOST_AUTO_TEST_CASE(metapop_test_sep_deme_and_ind_out_of_range)
{
    BOOST_REQUIRE_THROW(auto s = KTfwd::sample_separate(
                            pop, 2, std::vector<unsigned>({ 1000 }), true),
                        std::out_of_range);
}

BOOST_AUTO_TEST_CASE(metapop_test_deme_out_of_range)
{
    BOOST_REQUIRE_THROW(auto s = KTfwd::sample(rng.get(), pop, 2, 10, true),
                        std::out_of_range);
}

BOOST_AUTO_TEST_CASE(metapop_test_empty)
{
    auto s = KTfwd::sample(pop, 0, std::vector<unsigned>(), true);
    BOOST_REQUIRE(s.empty() == true);
}

BOOST_AUTO_TEST_CASE(metapop_test_ind_out_of_range)
{
    BOOST_REQUIRE_THROW(
        auto s = KTfwd::sample(pop, 2, std::vector<unsigned>({ 1000 }), true),
        std::out_of_range);
}

// Now, a metapop
BOOST_AUTO_TEST_CASE(metapop_1)
{
    pop.fixations.emplace_back(-0.1, 0, 0, 0);

    // Make sure that fiations end up in samples
    auto x = KTfwd::sample_separate(pop, 0, { 0, 1, 2 }, false);
    BOOST_REQUIRE_EQUAL(x.first.size(), 1);

    // Add a selected fixation
    pop.fixations.emplace_back(-0.2, 1.0, 0, 0);
    x = KTfwd::sample_separate(pop, 0, { 0, 1, 2 }, false);
    BOOST_REQUIRE_EQUAL(x.first.size(), 1);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(test_metapop_custom_sampling,
                         metapop_popgenmut_custom_fixture)

BOOST_AUTO_TEST_CASE(metapop_test_sep_deme_out_of_range_1)
{
    BOOST_REQUIRE_THROW(auto s
                        = KTfwd::sample_separate(rng.get(), pop, 2, 10, true),
                        std::out_of_range);
}

BOOST_AUTO_TEST_CASE(metapop_test_sep_deme_out_of_range_2)
{
    BOOST_REQUIRE_THROW(
        auto s = KTfwd::sample_separate(
            pop, 2, std::vector<unsigned>({ 10, 20, 50 }), true),
        std::out_of_range);
}

BOOST_AUTO_TEST_CASE(metapop_test_sep_throw)
{
    BOOST_REQUIRE_THROW(
        auto s = KTfwd::sample_separate(
            pop, 0, std::vector<unsigned>({ 10, 20, 50, 1000 }), true),
        std::out_of_range);
}

BOOST_AUTO_TEST_CASE(metapop_test_sep_empty)
{
    auto s = KTfwd::sample_separate(pop, 0, std::vector<unsigned>(), true);
    BOOST_REQUIRE(s.first.empty() == true);
    BOOST_REQUIRE(s.second.empty() == true);
}

BOOST_AUTO_TEST_CASE(metapop_test_sep_deme_and_ind_out_of_range)
{
    BOOST_REQUIRE_THROW(auto s = KTfwd::sample_separate(
                            pop, 2, std::vector<unsigned>({ 1000 }), true),
                        std::out_of_range);
}

BOOST_AUTO_TEST_CASE(metapop_test_deme_out_of_range)
{
    BOOST_REQUIRE_THROW(auto s = KTfwd::sample(rng.get(), pop, 2, 10, true),
                        std::out_of_range);
}

BOOST_AUTO_TEST_CASE(metapop_test_empty)
{
    auto s = KTfwd::sample(pop, 0, std::vector<unsigned>(), true);
    BOOST_REQUIRE(s.empty() == true);
}

BOOST_AUTO_TEST_CASE(metapop_test_ind_out_of_range)
{
    BOOST_REQUIRE_THROW(
        auto s = KTfwd::sample(pop, 2, std::vector<unsigned>({ 1000 }), true),
        std::out_of_range);
}

// Now, a metapop
BOOST_AUTO_TEST_CASE(metapop_1)
{
    pop.fixations.emplace_back(-0.1, 0, 0, 0);

    // Make sure that fiations end up in samples
    auto x = KTfwd::sample_separate(pop, 0, { 0, 1, 2 }, false);
    BOOST_REQUIRE_EQUAL(x.first.size(), 1);

    // Add a selected fixation
    pop.fixations.emplace_back(-0.2, 1.0, 0, 0);
    x = KTfwd::sample_separate(pop, 0, { 0, 1, 2 }, false);
    BOOST_REQUIRE_EQUAL(x.first.size(), 1);
}

BOOST_AUTO_TEST_SUITE_END()
