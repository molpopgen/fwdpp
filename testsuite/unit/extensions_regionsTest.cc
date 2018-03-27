/*
  \file extensions.cc
  API checks on fwdpp's extensions sub-library.
*/

#include <config.h>
#include <boost/test/unit_test.hpp>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/extensions/regions.hpp>
#include <fwdpp/type_traits.hpp>
#include <limits>
#include "../fixtures/sugar_fixtures.hpp"

using namespace fwdpp;

using poptype = slocuspop_popgenmut_fixture::poptype;
BOOST_FIXTURE_TEST_SUITE(test_extensions, slocuspop_popgenmut_fixture)

BOOST_AUTO_TEST_SUITE_END()

