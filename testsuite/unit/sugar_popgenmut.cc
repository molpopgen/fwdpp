// Test construction, etc., of mutation types
// in fwdpp/sugar

#include <sstream>
#include <boost/test/unit_test.hpp>
#include <fwdpp/types/mutation.hpp>

struct mutation_fixture
{
    fwdpp::mutation m;
    mutation_fixture() : m(0.1, -0.1, 0.5, 3, 1)
    {
    }
};

BOOST_AUTO_TEST_SUITE(test_mutation)

BOOST_FIXTURE_TEST_CASE(test_mutation_from_tuple, mutation_fixture)
{
    BOOST_REQUIRE_EQUAL(m.pos, 0.1);
    BOOST_REQUIRE_EQUAL(m.s, -0.1);
    BOOST_REQUIRE_EQUAL(m.h, 0.5);
    BOOST_REQUIRE_EQUAL(m.g, 3);
    BOOST_REQUIRE_EQUAL(m.xtra, 1);
    BOOST_REQUIRE_EQUAL(m.neutral, false);
}

BOOST_FIXTURE_TEST_CASE(test_mutation_comparison, mutation_fixture)
{
    auto m2 = m;
    BOOST_REQUIRE(m == m2);
    m.xtra = 17;
    BOOST_REQUIRE(!(m == m2));
}

BOOST_FIXTURE_TEST_CASE(test_serialize_mutation, mutation_fixture)
{
    std::ostringstream o;
    fwdpp::io::serialize_mutation<fwdpp::mutation>()(o, m);
    std::istringstream i(o.str());
    auto m2 = fwdpp::io::deserialize_mutation<fwdpp::mutation>()(i);

    BOOST_REQUIRE_EQUAL(m2.pos, 0.1);
    BOOST_REQUIRE_EQUAL(m2.s, -0.1);
    BOOST_REQUIRE_EQUAL(m2.h, 0.5);
    BOOST_REQUIRE_EQUAL(m2.g, 3);
    BOOST_REQUIRE_EQUAL(m2.xtra, 1);
    BOOST_REQUIRE_EQUAL(m2.neutral, false);
}

BOOST_AUTO_TEST_SUITE_END()
