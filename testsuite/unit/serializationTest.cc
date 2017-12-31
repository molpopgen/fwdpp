/*! \file serializationTest.cc
  \ingroup unit
  \brief Tests of low-level serialization functionality.
*/

/*
 * TODO: check how to test for error when writing to closed gzFile
 */

#include <config.h>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include "../fixtures/fwdpp_fixtures.hpp"
#include <fwdpp/io/scalar_serialization.hpp>
#include <boost/test/unit_test.hpp>

using gtype = fwdpp::gamete;

BOOST_AUTO_TEST_SUITE(test_scalar_writer)


BOOST_AUTO_TEST_SUITE_END()
