#include <iostream>
#include <fwdpp/GSLrng_t.hpp>
#include <boost/test/unit_test.hpp>
#include <gsl/gsl_randist.h>
#include <sstream>

BOOST_AUTO_TEST_SUITE(test_GSLrng)

BOOST_AUTO_TEST_CASE(is_moveable)
{
    fwdpp::GSLrng_t<fwdpp::GSL_RNG_MT19937> rng(101);

    auto rng2(std::move(rng));

    BOOST_CHECK(rng.get() == nullptr);
}

BOOST_AUTO_TEST_SUITE_END()
// BOOST_AUTO_TEST_CASE( serialize )
// {
//   fwdpp::GSLrng_t<fwdpp::GSL_RNG_MT19937> rng(101);

//   std::stringstream buffer;
//   rng.serialize(buffer);

//   fwdpp::GSLrng_t<fwdpp::GSL_RNG_MT19937> rng2(0); //make a new object with
//   different seed
//   //If all is cool, rng2 will have same state as rng 2
//   rng2.deserialize(buffer);
//   BOOST_CHECK_EQUAL(gsl_rng_get(rng.get()),gsl_rng_get(rng2.get()));
// }
