/*! \file demography.cc
  \ingroup unit
  \brief Testing functions in fwdpp/demography.hpp
*/

#include <config.h>
#include <algorithm>
#include <boost/test/unit_test.hpp>
#include "../fixtures/fwdpp_fixtures.hpp"
#include <fwdpp/debug.hpp>
#include <fwdpp/demography.hpp>

BOOST_FIXTURE_TEST_SUITE(test_demography, standard_empty_metapop_fixture)

BOOST_AUTO_TEST_CASE(low_level_copy_deme_test)
{
    diploids.emplace_back(dipvector_t(1000, std::make_pair(0, 0)));
    gametes.emplace_back(KTfwd::gamete(2000));
    BOOST_REQUIRE(KTfwd::check_sum(gametes, 2000));
    auto rv = KTfwd::copy_deme(mutations, mcounts, gametes, diploids,
                               0); // append a copy of deme 0 to the metapop
    BOOST_REQUIRE_EQUAL(rv, 0);
    BOOST_REQUIRE(KTfwd::check_sum(gametes, 4000));
    BOOST_REQUIRE_EQUAL(diploids.size(), 2);
    BOOST_REQUIRE(diploids[0] == diploids[1]);
    for (const auto &dips : diploids)
        {
            BOOST_REQUIRE(
                KTfwd::popdata_sane(dips, gametes, mutations, mcounts));
        }
}

BOOST_AUTO_TEST_CASE(low_level_copy_deme_test_out_of_range)
{
    diploids.emplace_back(dipvector_t(1000, std::make_pair(0, 0)));
    gametes.emplace_back(KTfwd::gamete(2000));
    auto rv
        = KTfwd::copy_deme(mutations, mcounts, gametes, diploids,
                           1); // Now, 1 is out of range, so we'll get an error
    BOOST_REQUIRE(rv != 0);
}

BOOST_AUTO_TEST_CASE(low_level_merge_deme_test)
{
    // Two-deme case
    {
        diploids = { dipvector_t(1000, std::make_pair(0, 0)),
                     dipvector_t(500, std::make_pair(1, 1)) };
        gametes = { KTfwd::gamete(2000), KTfwd::gamete(1000) };
        auto rv = KTfwd::merge_demes(diploids, 0, 1);
        BOOST_REQUIRE_EQUAL(rv, 0);
        BOOST_REQUIRE_EQUAL(diploids.size(), 1);
        BOOST_REQUIRE_EQUAL(diploids[0].size(), 1500);
        BOOST_REQUIRE(KTfwd::check_sum(gametes, 3000));
    }

    // Three-deme case
    {
        diploids = { dipvector_t(1000, std::make_pair(0, 0)),
                     dipvector_t(500, std::make_pair(1, 1)),
                     dipvector_t(250, std::make_pair(2, 2)) };
        gametes
            = { KTfwd::gamete(2000), KTfwd::gamete(1000), KTfwd::gamete(500) };
        auto rv = KTfwd::merge_demes(diploids, 1, 2);
        BOOST_REQUIRE_EQUAL(rv, 0);
        BOOST_REQUIRE_EQUAL(diploids.size(), 2);
        BOOST_REQUIRE_EQUAL(diploids[1].size(), 750);
        BOOST_REQUIRE(KTfwd::check_sum(gametes, 3500));
    }

    // Will give same result as previous test because 2 and 1 will get swapped
    {
        diploids = { dipvector_t(1000, std::make_pair(0, 0)),
                     dipvector_t(500, std::make_pair(1, 1)),
                     dipvector_t(250, std::make_pair(2, 2)) };
        gametes
            = { KTfwd::gamete(2000), KTfwd::gamete(1000), KTfwd::gamete(500) };
        auto rv = KTfwd::merge_demes(diploids, 2, 1);
        BOOST_REQUIRE_EQUAL(rv, 0);
        BOOST_REQUIRE_EQUAL(diploids.size(), 2);
        BOOST_REQUIRE_EQUAL(diploids[1].size(), 750);
        BOOST_REQUIRE(KTfwd::check_sum(gametes, 3500));
    }

    // Will not do anything b/c i == j
    {
        diploids = { dipvector_t(1000, std::make_pair(0, 0)),
                     dipvector_t(500, std::make_pair(1, 1)),
                     dipvector_t(250, std::make_pair(2, 2)) };
        gametes
            = { KTfwd::gamete(2000), KTfwd::gamete(1000), KTfwd::gamete(500) };
        auto rv = KTfwd::merge_demes(diploids, 1, 1);
        BOOST_REQUIRE_EQUAL(rv, 1);
    }

    // Two-deme case with index out of range
    {
        diploids = { dipvector_t(1000, std::make_pair(0, 0)),
                     dipvector_t(500, std::make_pair(1, 1)) };
        gametes = { KTfwd::gamete(2000), KTfwd::gamete(1000) };
        auto rv = KTfwd::merge_demes(diploids, 0, 2); // 2 is out of range
        BOOST_REQUIRE_EQUAL(rv, -1);
    }
}

BOOST_AUTO_TEST_CASE(low_level_remove_deme_test)
{
    {
        diploids = { dipvector_t(1000, std::make_pair(0, 0)),
                     dipvector_t(2000, std::make_pair(1, 1)),
                     dipvector_t(3000, std::make_pair(2, 2)) };
        gametes = { KTfwd::gamete(2000), KTfwd::gamete(4000),
                    KTfwd::gamete(6000) };

        // remove the second deme:
        auto rv = KTfwd::remove_deme(mutations, mcounts, gametes, diploids, 1);
        BOOST_REQUIRE(rv == 0);
        BOOST_REQUIRE(diploids.size() == 2);
        BOOST_REQUIRE(diploids[0].size() == 1000);
        BOOST_REQUIRE(diploids[1].size() == 3000);
        BOOST_REQUIRE(KTfwd::check_sum(gametes, 8000));
    }

    {
        diploids = { dipvector_t(1000, std::make_pair(0, 0)),
                     dipvector_t(2000, std::make_pair(1, 1)),
                     dipvector_t(3000, std::make_pair(2, 2)) };
        gametes = { KTfwd::gamete(2000), KTfwd::gamete(4000),
                    KTfwd::gamete(6000) };
        // deme index out of range:
        auto rv = KTfwd::remove_deme(mutations, mcounts, gametes, diploids, 3);
        BOOST_REQUIRE(rv == -1);
        BOOST_REQUIRE(KTfwd::check_sum(gametes, 12000));
    }
}

BOOST_AUTO_TEST_CASE(low_level_swap_demes_test)
{
    {
        std::size_t i = 0;
        diploids = { dipvector_t(1000, std::make_pair(i, i++)),
                     dipvector_t(2000, std::make_pair(i, i++)),
                     dipvector_t(3000, std::make_pair(i, i++)),
                     dipvector_t(4000, std::make_pair(i, i++)),
                     dipvector_t(5000, std::make_pair(i, i++)) };
        gametes
            = { KTfwd::gamete(2000), KTfwd::gamete(4000), KTfwd::gamete(6000),
                KTfwd::gamete(8000), KTfwd::gamete(10000) };
        auto rv = KTfwd::swap_demes(diploids, 0, 4);
        BOOST_REQUIRE_EQUAL(rv, 0);
        BOOST_REQUIRE(KTfwd::check_sum(gametes, 30000));
        BOOST_REQUIRE_EQUAL(diploids[0].size(), 5000);
        BOOST_REQUIRE_EQUAL(diploids[4].size(), 1000);
    }

    {
        std::size_t i = 0;
        diploids = { dipvector_t(1000, std::make_pair(i, i++)),
                     dipvector_t(2000, std::make_pair(i, i++)),
                     dipvector_t(3000, std::make_pair(i, i++)),
                     dipvector_t(4000, std::make_pair(i, i++)),
                     dipvector_t(5000, std::make_pair(i, i++)) };
        gametes
            = { KTfwd::gamete(2000), KTfwd::gamete(4000), KTfwd::gamete(6000),
                KTfwd::gamete(8000), KTfwd::gamete(10000) };
        auto rv = KTfwd::swap_demes(diploids, 1,
                                    1); // i==j, so nothing will happen
        BOOST_REQUIRE_EQUAL(rv, 1);
        BOOST_REQUIRE(KTfwd::check_sum(gametes, 30000));
    }

    {
        std::size_t i = 0;
        diploids = { dipvector_t(1000, std::make_pair(i, i++)),
                     dipvector_t(2000, std::make_pair(i, i++)),
                     dipvector_t(3000, std::make_pair(i, i++)),
                     dipvector_t(4000, std::make_pair(i, i++)),
                     dipvector_t(5000, std::make_pair(i, i++)) };
        gametes
            = { KTfwd::gamete(2000), KTfwd::gamete(4000), KTfwd::gamete(6000),
                KTfwd::gamete(8000), KTfwd::gamete(10000) };
        auto rv
            = KTfwd::swap_demes(diploids, 0, 5); // deme index j out of range
        BOOST_REQUIRE_EQUAL(rv, -1);
        BOOST_REQUIRE(KTfwd::check_sum(gametes, 30000));
    }
}

BOOST_AUTO_TEST_CASE(low_level_split_demes_test)
{
    gsl_rng *rng_ptr = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng_ptr, 0);

    {
        diploids = { dipvector_t(1000, std::make_pair(0, 0)) };
        gametes = { KTfwd::gamete(2000) };
        auto rv = KTfwd::split_deme(rng_ptr, mutations, mcounts, gametes,
                                    diploids, 0, 250, false);
        BOOST_REQUIRE_EQUAL(rv, 0);
        BOOST_REQUIRE_EQUAL(diploids.size(), 2);
        BOOST_REQUIRE_EQUAL(diploids[0].size(), 750);
        BOOST_REQUIRE_EQUAL(diploids[1].size(), 250);
        BOOST_REQUIRE(KTfwd::check_sum(gametes, 2000));
    }

    // with replacement
    {
        diploids = { dipvector_t(1000, std::make_pair(0, 0)) };
        gametes = { KTfwd::gamete(2000) };
        auto rv = KTfwd::split_deme(rng_ptr, mutations, mcounts, gametes,
                                    diploids, 0, 990, true);
        BOOST_REQUIRE_EQUAL(rv, 0);
        BOOST_REQUIRE_EQUAL(diploids.size(), 2);
        BOOST_REQUIRE_EQUAL(diploids[0].size(), 10);
        BOOST_REQUIRE_EQUAL(diploids[1].size(), 990);
        BOOST_REQUIRE(KTfwd::check_sum(gametes, 2000));
    }

    // index out of range
    {
        diploids = { dipvector_t(1000, std::make_pair(0, 0)) };
        gametes = { KTfwd::gamete(2000) };
        auto rv = KTfwd::split_deme(rng_ptr, mutations, mcounts, gametes,
                                    diploids, 1, 250, false);
        BOOST_REQUIRE_EQUAL(rv, -1);
    }
    gsl_rng_free(rng_ptr);
}

BOOST_AUTO_TEST_CASE(low_level_admix_demes_test)
{
    gsl_rng *rng_ptr = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng_ptr, 0);

    {
        diploids = { dipvector_t(1000, std::make_pair(0, 0)),
                     dipvector_t(1000, std::make_pair(1, 1)) };
        gametes = { KTfwd::gamete(2000), KTfwd::gamete(2000) };
        auto rv = KTfwd::admix_demes(rng_ptr, mutations, mcounts, gametes,
                                     diploids, 0, 1, 0.25, 1000, false);
        BOOST_REQUIRE_EQUAL(rv, 0);
        BOOST_REQUIRE(KTfwd::check_sum(gametes, 6000));
    }

    // with resampling
    {
        diploids = { dipvector_t(1000, std::make_pair(0, 0)),
                     dipvector_t(1000, std::make_pair(1, 1)) };
        gametes = { KTfwd::gamete(2000), KTfwd::gamete(2000) };
        auto rv = KTfwd::admix_demes(rng_ptr, mutations, mcounts, gametes,
                                     diploids, 0, 1, 0.25, 1000, true);
        BOOST_REQUIRE_EQUAL(rv, 0);
        BOOST_REQUIRE(KTfwd::check_sum(gametes, 6000));
    }

    // a pop index out of range
    {
        diploids = { dipvector_t(1000, std::make_pair(0, 0)),
                     dipvector_t(1000, std::make_pair(1, 1)) };
        gametes = { KTfwd::gamete(2000), KTfwd::gamete(2000) };
        auto rv = KTfwd::admix_demes(rng_ptr, mutations, mcounts, gametes,
                                     diploids, 0, 2, 0.25, 1000, false);
        BOOST_REQUIRE_EQUAL(rv, -1);
    }

    // Fraction of pop i ancestry < 0.
    {
        diploids = { dipvector_t(1000, std::make_pair(0, 0)),
                     dipvector_t(1000, std::make_pair(1, 1)) };
        gametes = { KTfwd::gamete(2000), KTfwd::gamete(2000) };
        auto rv = KTfwd::admix_demes(rng_ptr, mutations, mcounts, gametes,
                                     diploids, 0, 1, -0.25, 1000, false);
        BOOST_REQUIRE_EQUAL(rv, 1);
    }

    // Fraction of pop i ancestry >= 1.
    {
        diploids = { dipvector_t(1000, std::make_pair(0, 0)),
                     dipvector_t(1000, std::make_pair(1, 1)) };
        gametes = { KTfwd::gamete(2000), KTfwd::gamete(2000) };
        auto rv = KTfwd::admix_demes(rng_ptr, mutations, mcounts, gametes,
                                     diploids, 0, 1, 1., 1000, false);
        BOOST_REQUIRE_EQUAL(rv, 1);
    }

    /*
      When sampling w/o replacement, the new population size cannot be so large
      that the number of individuals required from deme i or j is >= than those
      demes'
      respective sizes.  This test samples so many individuals that this
      condition cannot be met,
      resulting in a return value of 1
    */

    {
        diploids = { dipvector_t(1000, std::make_pair(0, 0)),
                     dipvector_t(1000, std::make_pair(1, 1)) };
        gametes = { KTfwd::gamete(2000), KTfwd::gamete(2000) };
        // here, we want a new deme that is 50% ancestry from demes 0 and 1.
        // But, sampling w/o
        // replacement for 2000 individuals will require all individuals from
        // demes i and j.
        auto rv = KTfwd::admix_demes(rng_ptr, mutations, mcounts, gametes,
                                     diploids, 0, 1, 0.5, 2000, false);
        BOOST_REQUIRE_EQUAL(rv, 1);
    }

    // sampling w/replacement has no such limitations
    {
        diploids = { dipvector_t(1000, std::make_pair(0, 0)),
                     dipvector_t(1000, std::make_pair(1, 1)) };
        gametes = { KTfwd::gamete(2000), KTfwd::gamete(2000) };
        // here, we want a new deme that is 50% ancestry from demes 0 and 1.
        // But, sampling w/o
        // replacement for 2000 individuals will require all individuals from
        // demes i and j.
        auto rv = KTfwd::admix_demes(rng_ptr, mutations, mcounts, gametes,
                                     diploids, 0, 1, 0.5, 2000, true);
        BOOST_REQUIRE_EQUAL(rv, 0);
    }
    gsl_rng_free(rng_ptr);
}

// // #include <fwdpp/sugar/demography.hpp>

// // /*
// //   Now, we test the higher-level functions in the sugar sub-library.

// //   These functions update other data inside of KTfwd::metapop.

// //   These functions below are implemented via calls to the above functions.

// //   Thus, the tests below are API tests plus tests that
// KTfwd::sugar::metapop::Ns
// //   is getting properly updated..
// // */
// // BOOST_AUTO_TEST_CASE( copy_pop_test )
// // {
// //   {
// //     metapop_t mpop{1000};
// //     auto rv = KTfwd::copy_pop(mpop,0);
// //     BOOST_REQUIRE_EQUAL(rv,0);
// //     BOOST_REQUIRE_EQUAL(mpop.Ns.size(),2);
// //     BOOST_REQUIRE_EQUAL(mpop.Ns.size(),mpop.diploids.size());
// //     for(std::size_t i = 0 ; i < mpop.Ns.size() ; ++i)
// //       {
// // 	BOOST_REQUIRE_EQUAL(mpop.Ns[i],mpop.diploids[i].size());
// //       }
// //   }
// // }

// // BOOST_AUTO_TEST_CASE( merge_pops_test )
// // {
// //   {
// //     metapop_t mpop({1000,500});
// //     auto rv = KTfwd::merge_pops(mpop,0,1);
// //     BOOST_REQUIRE_EQUAL(rv,0);
// //     BOOST_REQUIRE(KTfwd::check_sum(mpop.gametes,3000));
// //     BOOST_REQUIRE_EQUAL(mpop.Ns.size(),1);
// //     BOOST_REQUIRE_EQUAL(mpop.Ns.size(),mpop.diploids.size());
// //     for(std::size_t i = 0 ; i < mpop.Ns.size() ; ++i)
// //       {
// // 	BOOST_REQUIRE_EQUAL(mpop.Ns[i],mpop.diploids[i].size());
// //       }
// //   }
// // }

// // BOOST_AUTO_TEST_CASE( swap_pops_test )
// // {
// //   {
// //     metapop_t mpop({1000,500});
// //     auto rv = KTfwd::swap_pops(mpop,0,1);
// //     BOOST_REQUIRE_EQUAL(rv,0);
// //     BOOST_REQUIRE(KTfwd::check_sum(mpop.gametes,3000));
// //     BOOST_REQUIRE_EQUAL(mpop.Ns.size(),2);
// //     BOOST_REQUIRE_EQUAL(mpop.Ns.size(),mpop.diploids.size());
// //     for(std::size_t i = 0 ; i < mpop.Ns.size() ; ++i)
// //       {
// // 	BOOST_REQUIRE_EQUAL(mpop.Ns[i],mpop.diploids[i].size());
// //       }
// //   }
// // }

// // BOOST_AUTO_TEST_CASE( remove_pop_test )
// // {
// //   {
// //     metapop_t mpop({1000,500});
// //     auto rv = KTfwd::remove_pop(mpop,1);
// //     BOOST_REQUIRE_EQUAL(rv,0);
// //     BOOST_REQUIRE(KTfwd::check_sum(mpop.gametes,2000));
// //     BOOST_REQUIRE_EQUAL(mpop.Ns.size(),1);
// //     BOOST_REQUIRE_EQUAL(mpop.Ns.size(),mpop.diploids.size());
// //     for(std::size_t i = 0 ; i < mpop.Ns.size() ; ++i)
// //       {
// // 	BOOST_REQUIRE_EQUAL(mpop.Ns[i],mpop.diploids[i].size());
// //       }
// //   }
// // }

// // BOOST_AUTO_TEST_CASE( split_pop_test )
// // {
// //   {
// //     metapop_t mpop{1000};
// //     KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);
// //     auto rv = KTfwd::split_pop(rng.get(),mpop,0,500,false);
// //     BOOST_REQUIRE_EQUAL(rv,0);
// //     BOOST_REQUIRE(KTfwd::check_sum(mpop.gametes,2000));
// //     BOOST_REQUIRE_EQUAL(mpop.Ns.size(),2);
// //     BOOST_REQUIRE_EQUAL(mpop.Ns.size(),mpop.diploids.size());
// //     for(std::size_t i = 0 ; i < mpop.Ns.size() ; ++i)
// //       {
// // 	BOOST_REQUIRE_EQUAL(mpop.Ns[i],mpop.diploids[i].size());
// //       }
// //   }
// //   {
// //     metapop_t mpop{1000};
// //     KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);
// //     auto rv = KTfwd::split_pop(rng.get(),mpop,0,500,true);
// //     BOOST_REQUIRE_EQUAL(rv,0);
// //     BOOST_REQUIRE(KTfwd::check_sum(mpop.gametes,2000));
// //     BOOST_REQUIRE_EQUAL(mpop.Ns.size(),2);
// //     BOOST_REQUIRE_EQUAL(mpop.Ns.size(),mpop.diploids.size());
// //     for(std::size_t i = 0 ; i < mpop.Ns.size() ; ++i)
// //       {
// // 	BOOST_REQUIRE_EQUAL(mpop.Ns[i],mpop.diploids[i].size());
// //       }
// //   }
// // }

// // BOOST_AUTO_TEST_CASE( admix_pops_test )
// // {
// //   {
// //     metapop_t mpop{1000,1000};
// //     KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);
// //     auto rv = KTfwd::admix_pops(rng.get(),mpop,0,1,0.5,1000);
// //     BOOST_REQUIRE_EQUAL(rv,0);
// //     BOOST_REQUIRE(KTfwd::check_sum(mpop.gametes,6000));
// //     BOOST_REQUIRE_EQUAL(mpop.Ns.size(),3);
// //     BOOST_REQUIRE_EQUAL(mpop.Ns.size(),mpop.diploids.size());
// //     for(std::size_t i = 0 ; i < mpop.Ns.size() ; ++i)
// //       {
// // 	BOOST_REQUIRE_EQUAL(mpop.Ns[i],mpop.diploids[i].size());
// //       }
// //   }
// // }

BOOST_AUTO_TEST_SUITE_END()
