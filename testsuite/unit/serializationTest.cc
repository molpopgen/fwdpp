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
#include <zlib.h>
#include "../fixtures/fwdpp_fixtures.hpp"
#include <fwdpp/internal/IOhelp.hpp>
#include <boost/test/unit_test.hpp>

using gtype = KTfwd::gamete;

BOOST_AUTO_TEST_SUITE(test_scalar_writer)

BOOST_AUTO_TEST_CASE(test_scalar_writer_exceptions)
/* The low-level function object KTfwd::fwdpp_internal::scalar_writer
 * will throw std::runtime_error from its operator() if an error is
 * encountered when writing.
 *
 * Note: attempting to write 0 bytes is considered an error as
 * far as this object is concerned!
 */
{
    std::ofstream foo;
    double x = 1.0;
    // Trying to write to invalid stream should throw
    BOOST_REQUIRE_THROW(KTfwd::fwdpp_internal::scalar_writer()(foo, &x),
                        std::runtime_error);

    std::ostringstream foo2;
    BOOST_REQUIRE_NO_THROW(KTfwd::fwdpp_internal::scalar_writer()(foo2, &x));
    // Should also throw when failbit, etc., is true
    foo2.setstate(std::ios::failbit);
    BOOST_REQUIRE_THROW(KTfwd::fwdpp_internal::scalar_writer()(foo2, &x),
                        std::runtime_error);

    gzFile gzf;
    // Will throw on an invalid gzfile
    BOOST_REQUIRE_THROW(KTfwd::fwdpp_internal::scalar_writer()(gzf, &x),
                        std::runtime_error);
    gzf = gzopen("test_scalar_writer_exceptions.gz", "wb");
    // Will throw if you attempt to write 0 bytes!!
    BOOST_REQUIRE_THROW(KTfwd::fwdpp_internal::scalar_writer()(gzf, &x, 0),
                        std::runtime_error);
    gzclose(gzf);
    // Throw when writing to closed file
    // BOOST_REQUIRE_THROW(KTfwd::fwdpp_internal::scalar_writer()(gzf,&x),std::runtime_error);
    // unlink("test_scalar_writer_exceptions.gz");
}

BOOST_AUTO_TEST_SUITE_END()
