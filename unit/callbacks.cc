#define BOOST_TEST_MODULE callbackTest
#define BOOST_TEST_DYN_LINK 

//This is really an API check

#include <config.h>
#include <fwdpp/extensions/callbacks.hpp>
#include <boost/test/unit_test.hpp>
#include <unistd.h>
#include <iterator>
#include <functional>
#include <list>
#include <iostream>

gsl_rng * r =  gsl_rng_alloc(gsl_rng_ranlxs2);

BOOST_AUTO_TEST_CASE(constant)
{
  gsl_rng_set(r,101);
  KTfwd::extensions::constant c(1);
}

BOOST_AUTO_TEST_CASE(exponential)
{
  gsl_rng_set(r,101);
  KTfwd::extensions::exponential e(1);
}

BOOST_AUTO_TEST_CASE(uniform)
{
  gsl_rng_set(r,101);
  KTfwd::extensions::uniform e(0,1);
}

BOOST_AUTO_TEST_CASE(beta)
{
  gsl_rng_set(r,101);
  KTfwd::extensions::beta b(1,2,1.0);
}

BOOST_AUTO_TEST_CASE(Gamma)
{
  gsl_rng_set(r,101);
  KTfwd::extensions::gamma b(1,2);
}

BOOST_AUTO_TEST_CASE(gaussian)
{
  gsl_rng_set(r,101);
  KTfwd::extensions::gaussian e(1);
}

BOOST_AUTO_TEST_CASE(shmodel)
{
  gsl_rng_set(r,101);
  //default construct
  KTfwd::extensions::shmodel x;
  //These both eval to the "right thing";
  x.s = std::bind(KTfwd::extensions::gamma(1,2),std::placeholders::_1);
  x.h = KTfwd::extensions::gamma(1,2);

  auto __s = x.s(r);
  auto __h = x.h(r);
}

BOOST_AUTO_TEST_CASE(shmodel2)
{
  gsl_rng_set(r,101);
  //Construct and consume the input arguments
  KTfwd::extensions::shmodel x(std::bind(KTfwd::extensions::gamma(1,2),std::placeholders::_1),
			       KTfwd::extensions::gamma(1,2));
  auto __s = x.s(r);
  auto __h = x.h(r);
}
