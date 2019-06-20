/*! \file serializationTest.cc
  \ingroup unit
  \brief Tests of low-level serialization functionality.
*/

#include <config.h>
#include <sstream>
#include <array>
#include <limits>
#include "../fixtures/fwdpp_fixtures.hpp"
#include <fwdpp/io/scalar_serialization.hpp>
#include <boost/test/unit_test.hpp>

using gtype = fwdpp::haploid_genome;

struct scalar_fixture
{
    std::array<double, 2> ad2, ad2_2;

    scalar_fixture()
        : ad2{ { -1, 1 } },
          ad2_2{ { std::numeric_limits<double>::quiet_NaN(),
                   std::numeric_limits<double>::quiet_NaN() } }
    {
    }
};

BOOST_AUTO_TEST_SUITE(test_scalar_writer)

BOOST_FIXTURE_TEST_CASE(serialize_std_array, scalar_fixture)
{
    fwdpp::io::scalar_writer writer;
    fwdpp::io::scalar_reader reader;
    std::stringstream buffer;
    writer(buffer,&ad2[0],2);
    reader(buffer,&ad2_2[0],2);
    BOOST_REQUIRE(ad2 == ad2_2);
}

BOOST_AUTO_TEST_SUITE_END()
