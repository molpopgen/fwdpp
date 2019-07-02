#include <config.h>
#include <cmath>
#include <fwdpp/genetic_map/poisson_interval.hpp>
#include <fwdpp/genetic_map/genetic_map.hpp>
#include <fwdpp/recbinder.hpp>
#include <boost/test/unit_test.hpp>
#include "../fixtures/rng_fixture.hpp"

void
apply_callback(const fwdpp::genetic_map_unit& gu, const gsl_rng* r)
{
    std::vector<double> breakpoints;
    gu(r, breakpoints);
}

BOOST_FIXTURE_TEST_SUITE(test_genetic_map, rng_fixture)

BOOST_AUTO_TEST_CASE(test_poisson_interval)
{
    BOOST_REQUIRE_NO_THROW(fwdpp::poisson_interval(0., 1., 1e-3));
    BOOST_REQUIRE_THROW(fwdpp::poisson_interval(0., 1., -1e-3),
                        std::invalid_argument);
    BOOST_REQUIRE_THROW(fwdpp::poisson_interval(1., 0., 1e-3),
                        std::invalid_argument);
    BOOST_REQUIRE_NO_THROW({
        fwdpp::poisson_interval p(0, 1, 1e-3);
        apply_callback(p, rng.get());
    });
}

BOOST_AUTO_TEST_CASE(test_genetic_map_moving_vector)
{
    BOOST_REQUIRE_NO_THROW({
        std::vector<std::unique_ptr<fwdpp::genetic_map_unit>> callbacks;
        callbacks.emplace_back(new fwdpp::poisson_interval(0., 1., 1e-3));
        fwdpp::genetic_map gm(std::move(callbacks));
        BOOST_REQUIRE_EQUAL(gm.size(), 1);
    });
}

BOOST_AUTO_TEST_CASE(test_genetic_map_binding)
{
    BOOST_REQUIRE_NO_THROW({
        std::vector<std::unique_ptr<fwdpp::genetic_map_unit>> callbacks;
        callbacks.emplace_back(new fwdpp::poisson_interval(0., 1., 1e-3));
        fwdpp::genetic_map gm(std::move(callbacks));
        auto r = fwdpp::recbinder(std::cref(gm), rng.get());
        auto breakpoints = r();
    });
}

BOOST_AUTO_TEST_SUITE_END()

