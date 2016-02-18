#define BOOST_TEST_MODULE GSLrng_t_test
#define BOOST_TEST_DYN_LINK 

#include <fwdpp/sugar/GSLrng_t.hpp>
#include <boost/test/unit_test.hpp>
#include <gsl/gsl_randist.h>

BOOST_AUTO_TEST_CASE( copy_construct_test )
{
  KTfwd::GSLrng_t<KTfwd::GSL_RNG_MT19937> rng(101);

  auto rng2(rng);

  //This tests that rng have same state
  BOOST_REQUIRE_EQUAL(gsl_rng_get(rng.get()),gsl_rng_get(rng2.get()));

  //If something went wrong with copy, destruction will fail here,
  //and we'll get a failure of the test
}
