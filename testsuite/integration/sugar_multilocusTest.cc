/*!
  \file sugar_multilocus.cc
  \ingroup unit
  \brief Testing KTfwd::multiloc
*/
#include <unistd.h>
#include <config.h>
#include <zlib.h>
#include <iostream>
#include <functional>
#include <algorithm>
#include <numeric>
#include <boost/test/unit_test.hpp>
#include <fwdpp/diploid.hh>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/sugar/multiloc.hpp>
#include <fwdpp/sugar/infsites.hpp>
#include <fwdpp/sugar/serialization.hpp>
#include "../fixtures/sugar_fixtures.hpp"
#include "../util/quick_evolve_sugar.hpp"

using poptype = multiloc_popgenmut_fixture::poptype;

void
simulate(multiloc_popgenmut_fixture &f, unsigned &generation)
{

    // Equal mutation and rec. rates per locus
    std::vector<double> mu(4, 0.005), rbw(3, 0.005);
    BOOST_REQUIRE_EQUAL(mu.size(), 4);
    for (; generation < 10; ++generation)
        {
            double wbar = KTfwd::sample_diploid(
                f.rng.get(), f.pop.gametes, f.pop.diploids, f.pop.mutations,
                f.pop.mcounts, 1000, &mu[0], f.mutmodels, f.recmodels, &rbw[0],
                [](const gsl_rng *__r, const double __d) {
                    return gsl_ran_binomial(__r, __d, 1);
                },
                std::bind(multiloc_popgenmut_fixture::multilocus_additive(),
                          std::placeholders::_1, std::placeholders::_2,
                          std::placeholders::_3),
                f.pop.neutral, f.pop.selected);
            assert(check_sum(f.pop.gametes, 8000));
            KTfwd::update_mutations(f.pop.mutations, f.pop.fixations,
                                    f.pop.fixation_times, f.pop.mut_lookup,
                                    f.pop.mcounts, generation, 2000);
        }
}

BOOST_AUTO_TEST_CASE(multiloc_sugar_test1)
{
    unsigned generation = 0;
    multiloc_popgenmut_fixture f(&generation);
    simulate(f, generation);
    poptype pop2(f.pop);
    BOOST_CHECK_EQUAL(f.pop == pop2, true);
}

BOOST_AUTO_TEST_CASE(multiloc_sugar_test2)
{
    unsigned generation = 0;
    multiloc_popgenmut_fixture f(&generation);
    simulate(f, generation);
    poptype pop2(0, 0);
    KTfwd::serialize s;
    std::stringstream buffer;
    s(buffer, f.pop, multiloc_popgenmut_fixture::mwriter());
    KTfwd::deserialize()(pop2, buffer, multiloc_popgenmut_fixture::mreader());
    BOOST_CHECK_EQUAL(f.pop == pop2, true);
}

BOOST_AUTO_TEST_CASE(multiloc_sugar_test2_gz)
{
    unsigned generation = 0;
    multiloc_popgenmut_fixture f(&generation);
    simulate(f, generation);
    poptype pop2(0, 0);
    gzFile gzf = gzopen("sugar_multilocus_out.gz", "wb");
    KTfwd::gzserialize()(gzf, f.pop, multiloc_popgenmut_fixture::mwriter());
    gzclose(gzf);
    gzf = gzopen("sugar_multilocus_out.gz", "rb");
    KTfwd::gzdeserialize()(pop2, gzf,
                           std::bind(multiloc_popgenmut_fixture::mreader(),
                                     std::placeholders::_1));
    gzclose(gzf);
    unlink("sugar_multilocus_out.gz");
    BOOST_CHECK_EQUAL(f.pop == pop2, true);
}

BOOST_AUTO_TEST_CASE(multiloc_sugar_test3)
{
    unsigned generation = 0;
    multiloc_popgenmut_fixture f(&generation);
    simulate(f, generation);
    poptype pop2(std::move(f.pop));
    BOOST_CHECK_EQUAL(f.pop == pop2, false);
}

BOOST_AUTO_TEST_CASE(multiloc_sugar_test4)
{
    unsigned generation = 0;
    multiloc_popgenmut_fixture f(&generation);
    simulate(f, generation);
    poptype pop2 = std::move(f.pop);
    BOOST_CHECK_EQUAL(f.pop == pop2, false);
}
