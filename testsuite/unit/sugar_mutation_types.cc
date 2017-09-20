// Test construction, etc., of mutation types
// in fwdpp/sugar

#include <boost/test/unit_test.hpp>
#include <fwdpp/sugar/popgenmut.hpp>
#include <fwdpp/sugar/generalmut.hpp>

struct popgenmut_tuple_wrapper
{
    KTfwd::popgenmut::constructor_tuple t;
    KTfwd::popgenmut m;
    popgenmut_tuple_wrapper() : t(std::make_tuple(0.1, -0.1, 0.5, 3, 1)), m(t)
    {
    }
};

BOOST_AUTO_TEST_SUITE(test_popgenmut)

    BOOST_FIXTURE_TEST_CASE(test_popgenmut_from_tuple, popgenmut_tuple_wrapper)
{
    BOOST_REQUIRE_EQUAL(m.pos, 0.1);
    BOOST_REQUIRE_EQUAL(m.s, -0.1);
    BOOST_REQUIRE_EQUAL(m.h, 0.5);
    BOOST_REQUIRE_EQUAL(m.g, 3);
    BOOST_REQUIRE_EQUAL(m.xtra, 1);
    BOOST_REQUIRE_EQUAL(m.neutral, false);
}

BOOST_AUTO_TEST_SUITE_END()
