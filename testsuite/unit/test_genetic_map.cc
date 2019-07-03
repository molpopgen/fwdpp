#include <config.h>
#include <cmath>
#include <iostream>
#include <fwdpp/genetic_map/poisson_interval.hpp>
#include <fwdpp/genetic_map/poisson_point.hpp>
#include <fwdpp/genetic_map/binomial_interval.hpp>
#include <fwdpp/genetic_map/binomial_point.hpp>
#include <fwdpp/genetic_map/fixed_number_crossovers.hpp>
#include <fwdpp/genetic_map/genetic_map.hpp>
#include <fwdpp/recbinder.hpp>
#include <boost/test/unit_test.hpp>
#include "../fixtures/rng_fixture.hpp"

std::vector<double>
apply_callback(const fwdpp::genetic_map_unit& gu, const gsl_rng* r)
{
    std::vector<double> breakpoints;
    gu(r, breakpoints);
    return breakpoints;
}

BOOST_FIXTURE_TEST_SUITE(test_poisson_interval, rng_fixture)

BOOST_AUTO_TEST_CASE(test_construction)
{
    BOOST_REQUIRE_NO_THROW(fwdpp::poisson_interval(0., 1., 1e-3));
    BOOST_REQUIRE_THROW(fwdpp::poisson_interval(0., 1., -1e-3),
                        std::invalid_argument);
    BOOST_REQUIRE_THROW(fwdpp::poisson_interval(1., 0., 1e-3),
                        std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(test_callback_nothrow)
{
    BOOST_REQUIRE_NO_THROW({
        fwdpp::poisson_interval p(0, 1, 1e-3);
        apply_callback(p, rng.get());
    });
}

BOOST_AUTO_TEST_CASE(test_clone_and_cast)
{
    BOOST_REQUIRE_NO_THROW({
        fwdpp::poisson_interval p(0, 1, 1e-3);
        auto c = p.clone();
        apply_callback(*c, rng.get());
    });
    fwdpp::poisson_interval p(0, 1, 1e-3);
    auto c = p.clone();
    BOOST_REQUIRE_EQUAL(c == nullptr, false);
    auto cast = dynamic_cast<fwdpp::poisson_interval*>(c.release());
    BOOST_REQUIRE_EQUAL(cast != nullptr, true);
    BOOST_REQUIRE_EQUAL(c == nullptr, true);
    BOOST_REQUIRE_EQUAL(cast->beg, p.beg);
    BOOST_REQUIRE_EQUAL(cast->end, p.end);
    BOOST_REQUIRE_EQUAL(cast->mean, p.mean);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(test_binomial_interval, rng_fixture)

BOOST_AUTO_TEST_CASE(test_construction)
{
    BOOST_REQUIRE_NO_THROW(fwdpp::binomial_interval(0., 1., 1e-3));
    BOOST_REQUIRE_THROW(fwdpp::binomial_interval(0., 1., -1e-3),
                        std::invalid_argument);
    BOOST_REQUIRE_THROW(fwdpp::binomial_interval(1., 0., 1e-3),
                        std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(test_callback_nothrow)
{
    BOOST_REQUIRE_NO_THROW({
        fwdpp::binomial_interval p(0, 1, 1);
        auto rv = apply_callback(p, rng.get());
        BOOST_REQUIRE_EQUAL(rv.size(), 1);
    });
}

BOOST_AUTO_TEST_CASE(test_clone_and_cast)
{
    BOOST_REQUIRE_NO_THROW({
        fwdpp::binomial_interval p(0, 1, 1e-3);
        auto c = p.clone();
        apply_callback(*c, rng.get());
    });
    fwdpp::binomial_interval p(0, 1, 1e-3);
    auto c = p.clone();
    BOOST_REQUIRE_EQUAL(c == nullptr, false);
    auto cast = dynamic_cast<fwdpp::binomial_interval*>(c.release());
    BOOST_REQUIRE_EQUAL(cast != nullptr, true);
    BOOST_REQUIRE_EQUAL(c == nullptr, true);
    BOOST_REQUIRE_EQUAL(cast->beg, p.beg);
    BOOST_REQUIRE_EQUAL(cast->end, p.end);
    BOOST_REQUIRE_EQUAL(cast->prob, p.prob);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(test_poisson_point, rng_fixture)

BOOST_AUTO_TEST_CASE(test_construction)
{
    BOOST_REQUIRE_NO_THROW(fwdpp::poisson_point(0., 1e-3));
    BOOST_REQUIRE_THROW(fwdpp::poisson_point(0., -1e-3),
                        std::invalid_argument);
    BOOST_REQUIRE_THROW(
        fwdpp::poisson_point(1., std::numeric_limits<double>::quiet_NaN()),
        std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(test_callback_nothrow)
{
    BOOST_REQUIRE_NO_THROW({
        fwdpp::poisson_point p(0, 1e-3);
        apply_callback(p, rng.get());
    });
}

BOOST_AUTO_TEST_CASE(test_clone_and_cast)
{
    BOOST_REQUIRE_NO_THROW({
        fwdpp::poisson_point p(0, 1e-3);
        auto c = p.clone();
        apply_callback(*c, rng.get());
    });
    fwdpp::poisson_point p(0, 1e-3);
    auto c = p.clone();
    BOOST_REQUIRE_EQUAL(c == nullptr, false);
    auto cast = dynamic_cast<fwdpp::poisson_point*>(c.release());
    BOOST_REQUIRE_EQUAL(cast != nullptr, true);
    BOOST_REQUIRE_EQUAL(c == nullptr, true);
    BOOST_REQUIRE_EQUAL(cast->position, p.position);
    BOOST_REQUIRE_EQUAL(cast->mean, p.mean);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(test_binomial_point, rng_fixture)

BOOST_AUTO_TEST_CASE(test_construction)
{
    BOOST_REQUIRE_NO_THROW(fwdpp::binomial_point(0., 1e-3));
    BOOST_REQUIRE_THROW(fwdpp::binomial_point(0., -1e-3),
                        std::invalid_argument);
    BOOST_REQUIRE_THROW(
        fwdpp::binomial_point(1., std::numeric_limits<double>::quiet_NaN()),
        std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(test_callback_nothrow)
{
    BOOST_REQUIRE_NO_THROW({
        fwdpp::binomial_point p(0, 1);
        auto rv = apply_callback(p, rng.get());
        BOOST_REQUIRE_EQUAL(rv.size(), 1);
        BOOST_REQUIRE_EQUAL(rv.front(), 0);
    });
}

BOOST_AUTO_TEST_CASE(test_clone_and_cast)
{
    BOOST_REQUIRE_NO_THROW({
        fwdpp::binomial_point p(0, 1e-3);
        auto c = p.clone();
        apply_callback(*c, rng.get());
    });
    fwdpp::binomial_point p(0, 1e-3);
    auto c = p.clone();
    BOOST_REQUIRE_EQUAL(c == nullptr, false);
    auto cast = dynamic_cast<fwdpp::binomial_point*>(c.release());
    BOOST_REQUIRE_EQUAL(cast != nullptr, true);
    BOOST_REQUIRE_EQUAL(c == nullptr, true);
    BOOST_REQUIRE_EQUAL(cast->position, p.position);
    BOOST_REQUIRE_EQUAL(cast->prob, p.prob);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(test_fixed_number_crossovers, rng_fixture)

BOOST_AUTO_TEST_CASE(test_construction)
{
    BOOST_REQUIRE_NO_THROW(fwdpp::fixed_number_crossovers(0, 1, 10));
    BOOST_REQUIRE_THROW(fwdpp::fixed_number_crossovers(1, 0, 10),
                        std::invalid_argument);
    BOOST_REQUIRE_THROW(fwdpp::fixed_number_crossovers(0, 1, -1),
                        std::invalid_argument);
    BOOST_REQUIRE_THROW(fwdpp::fixed_number_crossovers(
                            std::numeric_limits<double>::quiet_NaN(), 1, 1),
                        std::invalid_argument);
    BOOST_REQUIRE_THROW(fwdpp::fixed_number_crossovers(
                            1., std::numeric_limits<double>::quiet_NaN(), 1),
                        std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(test_callback_nothrow)
{
    BOOST_REQUIRE_NO_THROW({
        fwdpp::fixed_number_crossovers p(0, 1, 5);
        auto rv = apply_callback(p, rng.get());
        BOOST_REQUIRE_EQUAL(rv.size(), 5);
    });
}

BOOST_AUTO_TEST_CASE(test_clone_and_cast)
{
    BOOST_REQUIRE_NO_THROW({
        fwdpp::fixed_number_crossovers p(0, 1, 1);
        auto c = p.clone();
        apply_callback(*c, rng.get());
    });
    fwdpp::fixed_number_crossovers p(0, 1, 1);
    auto c = p.clone();
    BOOST_REQUIRE_EQUAL(c == nullptr, false);
    auto cast = dynamic_cast<fwdpp::fixed_number_crossovers*>(c.release());
    BOOST_REQUIRE_EQUAL(cast != nullptr, true);
    BOOST_REQUIRE_EQUAL(c == nullptr, true);
    BOOST_REQUIRE_EQUAL(cast->beg, p.beg);
    BOOST_REQUIRE_EQUAL(cast->end, p.end);
    BOOST_REQUIRE_EQUAL(cast->nxovers, p.nxovers);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(test_genetic_map, rng_fixture)

BOOST_AUTO_TEST_CASE(test_moving_vector)
{
    BOOST_REQUIRE_NO_THROW({
        std::vector<std::unique_ptr<fwdpp::genetic_map_unit>> callbacks;
        callbacks.emplace_back(new fwdpp::poisson_interval(0., 1., 1e-3));
        fwdpp::genetic_map gm(std::move(callbacks));
        BOOST_REQUIRE_EQUAL(gm.size(), 1);
        auto r = fwdpp::recbinder(std::cref(gm), rng.get());
        auto breakpoints = r();
    });
}

BOOST_AUTO_TEST_CASE(test_adding_poisson_interval)
{
    fwdpp::poisson_interval p(0, 1, 1e-1);
    fwdpp::genetic_map gm;
    gm.add_callback(p);
    BOOST_REQUIRE_EQUAL(gm.size(), 1);
    auto r = fwdpp::recbinder(std::cref(gm), rng.get());
    auto breakpoints = r();
}

BOOST_AUTO_TEST_CASE(test_adding_poisson_point)
{
    fwdpp::poisson_point p(0, 1);
    fwdpp::genetic_map gm;
    gm.add_callback(p);
    BOOST_REQUIRE_EQUAL(gm.size(), 1);
    auto r = fwdpp::recbinder(std::cref(gm), rng.get());
    auto breakpoints = r();
}

BOOST_AUTO_TEST_CASE(test_adding_binomial_point)
{
    fwdpp::binomial_point p(0, 1);
    fwdpp::genetic_map gm;
    gm.add_callback(p);
    BOOST_REQUIRE_EQUAL(gm.size(), 1);
    auto r = fwdpp::recbinder(std::cref(gm), rng.get());
    auto breakpoints = r();
}

BOOST_AUTO_TEST_SUITE_END()

