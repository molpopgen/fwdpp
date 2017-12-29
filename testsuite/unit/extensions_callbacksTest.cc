// This is really an API check

#include <config.h>
#include <cmath>
#include <fwdpp/extensions/callbacks.hpp>
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(test_extensions_callbacks)

gsl_rng *r = gsl_rng_alloc(gsl_rng_ranlxs2);

BOOST_AUTO_TEST_CASE(constant) { fwdpp::extensions::constant c(1); }

BOOST_AUTO_TEST_CASE(exponential) { fwdpp::extensions::exponential e(1); }

BOOST_AUTO_TEST_CASE(uniform) { fwdpp::extensions::uniform e(0, 1); }

BOOST_AUTO_TEST_CASE(beta) { fwdpp::extensions::beta b(1, 2, 1.0); }

BOOST_AUTO_TEST_CASE(Gamma) { fwdpp::extensions::gamma b(1, 2); }

BOOST_AUTO_TEST_CASE(gaussian) { fwdpp::extensions::gaussian e(1); }

BOOST_AUTO_TEST_CASE(shmodel)
{
    gsl_rng_set(r, 101);
    // default construct
    fwdpp::extensions::shmodel x;
    // These both eval to the "right thing";
    x.s = std::bind(fwdpp::extensions::gamma(1, 2), std::placeholders::_1);
    x.h = fwdpp::extensions::gamma(1, 2);

    auto __s = x.s(r);
    auto __h = x.h(r);
    BOOST_REQUIRE_EQUAL(std::isfinite(__s), true);
    BOOST_REQUIRE_EQUAL(std::isfinite(__h), true);
}

BOOST_AUTO_TEST_CASE(shmodel2)
{
    gsl_rng_set(r, 101);
    // Construct and consume the input arguments
    fwdpp::extensions::shmodel x(
        std::bind(fwdpp::extensions::gamma(1, 2), std::placeholders::_1),
        fwdpp::extensions::gamma(1, 2));
    auto __s = x.s(r);
    auto __h = x.h(r);
    BOOST_REQUIRE_EQUAL(std::isfinite(__s), true);
    BOOST_REQUIRE_EQUAL(std::isfinite(__h), true);
}

BOOST_AUTO_TEST_CASE(point_mass_test1)
{
    gsl_rng_set(r, 101);
    auto x = fwdpp::extensions::uniform(-1., -1.);
    auto y = x(r);
    BOOST_REQUIRE_EQUAL(y, -1.0);

    // Free here--kinda lame..
    gsl_rng_free(r);
}

// The callbacks can throw exceptions if their parameters aren't valid

BOOST_AUTO_TEST_CASE(callback_exceptions)
{
    {
        // inf
        BOOST_REQUIRE_THROW(fwdpp::extensions::constant(1. / 0.),
                            std::invalid_argument);
    }

    {
        // nan
        BOOST_REQUIRE_THROW(fwdpp::extensions::constant(std::nan("")),
                            std::invalid_argument);
    }

    {
        // first arg not finite
        BOOST_REQUIRE_THROW(fwdpp::extensions::uniform(1. / 0., 1.),
                            std::invalid_argument);
    }

    {
        // 2nd arg not finite
        BOOST_REQUIRE_THROW(fwdpp::extensions::uniform(1., 1. / 0.),
                            std::invalid_argument);
    }

    {
        // min > max
        BOOST_REQUIRE_THROW(fwdpp::extensions::uniform(1., 0.99),
                            std::invalid_argument);
    }

    {
        // a not finite
        BOOST_REQUIRE_THROW(fwdpp::extensions::beta(std::nan(""), 1.),
                            std::invalid_argument);
    }

    {
        // b not finite
        BOOST_REQUIRE_THROW(fwdpp::extensions::beta(1., std::nan("")),
                            std::invalid_argument);
    }

    {
        // f not finite
        BOOST_REQUIRE_THROW(fwdpp::extensions::beta(1., 1., std::nan("")),
                            std::invalid_argument);
    }

    {
        // a <= 0.
        BOOST_REQUIRE_THROW(fwdpp::extensions::beta(0., 1.),
                            std::invalid_argument);
    }

    {
        // b <= 0.
        BOOST_REQUIRE_THROW(fwdpp::extensions::beta(1., 0.),
                            std::invalid_argument);
    }

    {
        // f <= 0.
        BOOST_REQUIRE_THROW(fwdpp::extensions::beta(1., 1., 0.),
                            std::invalid_argument);
    }

    {
        // sd = 0
        BOOST_REQUIRE_THROW(fwdpp::extensions::gaussian(0.),
                            std::invalid_argument);
    }

    {
        // sd < 0
        BOOST_REQUIRE_THROW(fwdpp::extensions::gaussian(-1e-6),
                            std::invalid_argument);
    }

    {
        // sd not finite
        BOOST_REQUIRE_THROW(fwdpp::extensions::gaussian(1. / 0.),
                            std::invalid_argument);
    }

    {
        // mean not finite
        BOOST_REQUIRE_THROW(fwdpp::extensions::gamma(1. / 0., 1.),
                            std::invalid_argument);
    }

    {
        // shape not finite
        BOOST_REQUIRE_THROW(fwdpp::extensions::gamma(1.0, 1. / 0.),
                            std::invalid_argument);
    }

    {
        //!(shape>0)
        BOOST_REQUIRE_THROW(fwdpp::extensions::gamma(1.0, 0.),
                            std::invalid_argument);
    }
}

BOOST_AUTO_TEST_CASE(vector_shmodel)
{
    using namespace fwdpp;
    // PS, uniform initialization rocks...
    std::vector<extensions::shmodel> callbacks{
        { extensions::constant(1.), extensions::constant(0.) },
        { extensions::exponential(1.), extensions::exponential(1.) },
        { extensions::uniform(1., 2.), extensions::uniform(1., 2.) },
        { extensions::beta(1., 2.),
          extensions::beta(1., 2.) }, // defaults to factor = 1
        { extensions::beta(1., 2., 0.25),
          extensions::beta(1., 2., 0.25) }, // pass all 3 params to constructor
        { extensions::gaussian(1.), extensions::gaussian(1.) },
        { extensions::gamma(1., 0.1), extensions::gamma(1., 0.1) }
    };
}

BOOST_AUTO_TEST_SUITE_END()
