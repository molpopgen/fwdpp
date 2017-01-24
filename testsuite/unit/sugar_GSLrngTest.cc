#include <iostream>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <boost/test/unit_test.hpp>
#include <gsl/gsl_randist.h>
#include <sstream>

BOOST_AUTO_TEST_SUITE(test_GSLrng)

BOOST_AUTO_TEST_CASE(copy_construct_test)
{
    KTfwd::GSLrng_t<KTfwd::GSL_RNG_MT19937> rng(101);

    auto rng2(rng);

    // This tests that rng have same state
    BOOST_REQUIRE_EQUAL(gsl_rng_get(rng.get()), gsl_rng_get(rng2.get()));

    // If something went wrong with copy, destruction will fail here,
    // and we'll get a failure of the test
}

BOOST_AUTO_TEST_CASE(is_moveable)
{
    KTfwd::GSLrng_t<KTfwd::GSL_RNG_MT19937> rng(101);

    auto rng2(std::move(rng));

    BOOST_CHECK(rng.get() == nullptr);
}

BOOST_AUTO_TEST_SUITE_END()
// BOOST_AUTO_TEST_CASE( serialize )
// {
//   KTfwd::GSLrng_t<KTfwd::GSL_RNG_MT19937> rng(101);

//   std::stringstream buffer;
//   rng.serialize(buffer);

//   KTfwd::GSLrng_t<KTfwd::GSL_RNG_MT19937> rng2(0); //make a new object with
//   different seed
//   //If all is cool, rng2 will have same state as rng 2
//   rng2.deserialize(buffer);
//   BOOST_CHECK_EQUAL(gsl_rng_get(rng.get()),gsl_rng_get(rng2.get()));
// }
