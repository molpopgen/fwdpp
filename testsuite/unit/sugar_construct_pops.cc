// Test the construction of populations
// from user input data

#include <iostream>
#include <boost/test/unit_test.hpp>
#include <fwdpp/sugar/singlepop.hpp>
#include <fwdpp/sugar/popgenmut.hpp>
#include <testsuite/fixtures/sugar_fixtures.hpp>

struct singlepop_objects
{
    using poptype = singlepop_popgenmut_fixture::poptype;

    poptype::dipvector_t diploids;
    poptype::mcont_t mutations;
    poptype::gcont_t gametes;

    using mutation_container
        = poptype::gcont_t::value_type::mutation_container;

    singlepop_objects() : diploids{}, mutations{}, gametes{}
    {
        // Add some mutations
        for (unsigned i = 0; i < 3; ++i)
            {
                // position = i, effect size = i
                // meaning muts 1 and 2 not neutral
                mutations.emplace_back(i, i, 1, 0);
            }

        // Add two gametes
        gametes.emplace_back(1, mutation_container{ 0 },
                             mutation_container{ 1 });
        gametes.emplace_back(1, mutation_container{}, mutation_container{ 2 });

        // Add a diploid
        diploids.emplace_back(0, 1);
    }
};

BOOST_AUTO_TEST_SUITE(test_construct_pops)

BOOST_FIXTURE_TEST_CASE(test_singlepop, singlepop_objects)
{
    BOOST_REQUIRE_EQUAL(mutations[0].neutral, true);
    BOOST_REQUIRE_EQUAL(mutations[1].neutral, false);
    BOOST_REQUIRE_EQUAL(mutations[2].neutral, false);

    BOOST_REQUIRE_NO_THROW(poptype p(diploids, gametes, mutations););
}

BOOST_FIXTURE_TEST_CASE(test_singlepop_bad_gamete_counts, singlepop_objects)
{
    // Now, make incorrect data to trigger exceptions
    diploids[0].second = 0;
    BOOST_REQUIRE_THROW(poptype p(diploids, gametes, mutations),
                        std::runtime_error);

    // restore
    diploids[0].second = 1;
    // make the gamete data incorrect
    gametes[0].n++;
    BOOST_REQUIRE_THROW(poptype p(diploids, gametes, mutations),
                        std::runtime_error);
}

BOOST_FIXTURE_TEST_CASE(test_singlepop_incorrect_mutation_storage, singlepop_objects)
{
    std::swap(gametes[0].mutations,gametes[0].smutations);
    BOOST_REQUIRE_THROW(poptype p(diploids, gametes, mutations),
                        std::logic_error);
}

BOOST_AUTO_TEST_SUITE_END()
