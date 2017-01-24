/*! \file sugar_metapop.cc
  \ingroup unit
  \brief Testing KTfwd::metapop
*/
#include <config.h>
#include <boost/test/unit_test.hpp>
#include "../fixtures/sugar_fixtures.hpp"
#include "../util/quick_evolve_sugar.hpp"

using poptype = metapop_popgenmut_fixture::poptype;
using spoptype = singlepop_popgenmut_fixture::poptype;

BOOST_AUTO_TEST_CASE(uniform_init)
// Make sure C++11-style initialization is working ok
{
    {
        poptype mpop{ 100, 500 };
        BOOST_REQUIRE(mpop.diploids.size() == 2);
        BOOST_REQUIRE(mpop.diploids[0].size() == 100);
        BOOST_REQUIRE(mpop.diploids[1].size() == 500);
    }

    {
        poptype mpop({ 100, 500 });
        BOOST_REQUIRE(mpop.diploids.size() == 2);
        BOOST_REQUIRE(mpop.diploids[0].size() == 100);
        BOOST_REQUIRE(mpop.diploids[1].size() == 500);
    }
}

BOOST_AUTO_TEST_CASE(copy_construct_metapop_from_singlepop)
{
    spoptype pop(1000);
    BOOST_REQUIRE_EQUAL(pop.N, 1000);
    poptype mpop(pop);
    BOOST_REQUIRE_EQUAL(mpop.Ns.size(), 1);
    BOOST_REQUIRE_EQUAL(mpop.Ns[0], 1000);
    BOOST_REQUIRE_EQUAL(mpop.diploids.size(), 1);
    BOOST_REQUIRE_EQUAL(mpop.diploids[0].size(), 1000);
    BOOST_REQUIRE_EQUAL(pop.mutations == mpop.mutations, true);
    BOOST_REQUIRE_EQUAL(pop.gametes == mpop.gametes, true);
    BOOST_REQUIRE_EQUAL(pop.diploids == mpop.diploids[0], true);
}

BOOST_AUTO_TEST_CASE(move_construct_metapop_from_singlepop)
{
    spoptype pop(1000);
    BOOST_REQUIRE_EQUAL(pop.N, 1000);
    poptype mpop(std::move(pop));
    BOOST_REQUIRE_EQUAL(mpop.Ns.size(), 1);
    BOOST_REQUIRE_EQUAL(mpop.Ns[0], 1000);
    BOOST_REQUIRE_EQUAL(mpop.diploids.size(), 1);
    BOOST_REQUIRE_EQUAL(mpop.diploids[0].size(), 1000);

    // We do not test sizes of elements in pop b/c
    // what happens will be a little bit compiler-dependent.
    // Typically, though, most or all of the elements in pop
    // should be empty.
}

BOOST_FIXTURE_TEST_SUITE(test_metapop_copy_move, metapop_popgenmut_fixture)

BOOST_AUTO_TEST_CASE(metapop_sugar_test1)
{
    simulate_metapop(pop);
    auto pop2(pop);
    BOOST_CHECK_EQUAL(pop == pop2, true);
}

BOOST_AUTO_TEST_CASE(metapop_sugar_test3)
{
    simulate_metapop(pop);
    auto pop2(std::move(pop));
    BOOST_CHECK_EQUAL(pop == pop2, false);
}

BOOST_AUTO_TEST_CASE(metapop_sugar_test4)
{
    simulate_metapop(pop);
    auto pop2 = std::move(pop);
    BOOST_CHECK_EQUAL(pop == pop2, false);
}

BOOST_AUTO_TEST_SUITE_END()

/*
  These next two derived classes mimic what software
  packages like fwdpy do, which is to extend sugar types
  with data that they need.
*/

struct spop_derived : public spoptype
{
    unsigned generation;
    spop_derived(unsigned N) : spoptype(N), generation(0) {}
};

struct mpop_derived : public poptype
{
    unsigned generation;
    mpop_derived(std::initializer_list<unsigned> &Ns)
        : poptype(Ns), generation(0)
    {
    }
    mpop_derived(const spop_derived &s) : poptype(s), generation(s.generation)
    {
    }
    mpop_derived(spop_derived &&s)
        : poptype(std::move(s)), generation(s.generation)
    {
    }
};

BOOST_AUTO_TEST_SUITE(test_derive_from_sugar_types)

BOOST_AUTO_TEST_CASE(copy_construct_mpop_derived_from_spop_derived)
{
    spop_derived pop(1000);
    BOOST_REQUIRE_EQUAL(pop.N, 1000);
    mpop_derived mpop(pop);
    BOOST_REQUIRE_EQUAL(mpop.Ns.size(), 1);
    BOOST_REQUIRE_EQUAL(mpop.Ns[0], 1000);
    BOOST_REQUIRE_EQUAL(mpop.diploids.size(), 1);
    BOOST_REQUIRE_EQUAL(mpop.diploids[0].size(), 1000);
    BOOST_REQUIRE_EQUAL(pop.mutations == mpop.mutations, true);
    BOOST_REQUIRE_EQUAL(pop.gametes == mpop.gametes, true);
    BOOST_REQUIRE_EQUAL(pop.diploids == mpop.diploids[0], true);
}

BOOST_AUTO_TEST_SUITE_END()
