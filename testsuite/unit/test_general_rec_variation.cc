/*!
  \file test_general_rec_variation.cc
  \ingroup unit
  \brief Unit tests for KTfwd::general_rec_variation.  Also tests the built-in policies accompanying this type.
*/

#include <boost/test/unit_test.hpp>
#include <fwdpp/general_rec_variation.hpp>
#include <config.h>

struct fixture
{
    const gsl_rng* r;
    const double not_a_num, inf;
    fixture()
        : r{ gsl_rng_alloc(gsl_rng_taus) },
          not_a_num{ std::numeric_limits<double>::quiet_NaN() },
          inf{ std::numeric_limits<double>::infinity() }
    {
        gsl_rng_set(r, 42);
    }
};

BOOST_FIXTURE_TEST_SUITE(test_general_rec_variation, fixture)

BOOST_AUTO_TEST_CASE(test_poisson_interval)
{
    BOOST_REQUIRE_THROW(KTfwd::poisson_interval(NULL, 1e-3, 0., 1.),
                        std::invalid_argument);
    BOOST_REQUIRE_THROW(KTfwd::poisson_interval(r, not_a_num, 0., 1.),
                        std::invalid_argument);
    BOOST_REQUIRE_THROW(KTfwd::poisson_interval(r, 1e-3, not_a_num, 1.),
                        std::invalid_argument);
    BOOST_REQUIRE_THROW(KTfwd::poisson_interval(r, 1e-3, 0., not_a_num),
                        std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(test_crossover_point)
{
    BOOST_REQUIRE_THROW(KTfwd::crossover_point(NULL, 1e-3, -2),
                        std::invalid_argument);
    BOOST_REQUIRE_THROW(KTfwd::crossover_point(r, inf, -2),
                        std::invalid_argument);
    BOOST_REQUIRE_THROW(KTfwd::crossover_point(r, -2, not_a_num),
                        std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(test_basic_use)
{
    KTfwd::general_rec_variation rv;

    rv.recmap.push_back(KTfwd::poisson_interval(r, 1e-3, 0., 1.));
    rv.recmap.push_back(KTfwd::crossover_point(r, 1e-3, 0.5));
    rv.recmap.push_back(KTfwd::crossover_point(r, 1e-3, 0.5, true));
    rv.recmap.push_back(KTfwd::crossover_point(r, 1e-3, 0.5, false));

    auto x = rv();
}

BOOST_AUTO_TEST_SUITE_END()
