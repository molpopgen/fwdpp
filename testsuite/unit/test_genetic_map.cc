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

BOOST_AUTO_TEST_SUITE(test_cast_function)

BOOST_AUTO_TEST_CASE(test_continuous)
{
    auto f = fwdpp::detail::generate_genetic_map_unit_cast_function(false);
    auto x = f(10.2);
    BOOST_REQUIRE_EQUAL(x, 10.2);
}

BOOST_AUTO_TEST_CASE(test_discrete)
{
    auto f = fwdpp::detail::generate_genetic_map_unit_cast_function(true);
    auto x = f(10.2);
    BOOST_REQUIRE_EQUAL(x, 10.0);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(test_poisson_interval, rng_fixture)

BOOST_AUTO_TEST_CASE(test_continuous_construction)
{
    BOOST_REQUIRE_NO_THROW(fwdpp::poisson_interval(0., 1., 1e-3, false));
    BOOST_REQUIRE_THROW(fwdpp::poisson_interval(0., 1., -1e-3, false),
                        std::invalid_argument);
    BOOST_REQUIRE_THROW(fwdpp::poisson_interval(1., 0., 1e-3, false),
                        std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(test_continuous_init)
{
    fwdpp::poisson_interval p(0.0, 10.0, 5.0, false);
    BOOST_REQUIRE(p.discrete() == false);
}

BOOST_AUTO_TEST_CASE(test_callback_nothrow)
{
    BOOST_REQUIRE_NO_THROW({
        fwdpp::poisson_interval p(0, 1, 1e-3, false);
        apply_callback(p, rng.get());
    });
}

BOOST_AUTO_TEST_CASE(test_clone_and_cast)
{
    BOOST_REQUIRE_NO_THROW({
        fwdpp::poisson_interval p(0, 1, 1e-3, false);
        auto c = p.clone();
        apply_callback(*c, rng.get());
    });
    fwdpp::poisson_interval p(0, 1, 1e-3, false);
    auto c = p.clone();
    BOOST_REQUIRE_EQUAL(c == nullptr, false);
    auto cast = dynamic_cast<fwdpp::poisson_interval*>(c.get());
    BOOST_REQUIRE_EQUAL(cast != nullptr, true);
    BOOST_REQUIRE_EQUAL(cast->beg, p.beg);
    BOOST_REQUIRE_EQUAL(cast->end, p.end);
    BOOST_REQUIRE_EQUAL(cast->mean, p.mean);
    BOOST_REQUIRE_EQUAL(cast->discrete(), p.discrete());
}

BOOST_AUTO_TEST_CASE(test_discrete_construction)
{
    fwdpp::poisson_interval p(0.0, 10.0, 5.0, true);
    BOOST_REQUIRE(p.discrete());
}

BOOST_FIXTURE_TEST_CASE(test_discrete_call_operator, rng_fixture)
{
    fwdpp::poisson_interval p(0.0, 10.0, 5.0, true);
    auto rv = apply_callback(p, rng.get());
    BOOST_REQUIRE(!rv.empty());
    for (auto i : rv)
        {
            BOOST_REQUIRE_EQUAL(i, std::floor(i));
        }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(test_binomial_interval, rng_fixture)

BOOST_AUTO_TEST_CASE(test_continuous_construction)
{
    BOOST_REQUIRE_NO_THROW(fwdpp::binomial_interval(0., 1., 1e-3, false));
    BOOST_REQUIRE_THROW(fwdpp::binomial_interval(0., 1., -1e-3, false),
                        std::invalid_argument);
    BOOST_REQUIRE_THROW(fwdpp::binomial_interval(1., 0., 1e-3, false),
                        std::invalid_argument);
    fwdpp::binomial_interval b(0., 1., 1e-3, false);
    BOOST_REQUIRE_EQUAL(b.discrete(), false);
}

BOOST_AUTO_TEST_CASE(test_callback_nothrow)
{
    BOOST_REQUIRE_NO_THROW({
        fwdpp::binomial_interval p(0, 1, 1, false);
        auto rv = apply_callback(p, rng.get());
        BOOST_REQUIRE_EQUAL(rv.size(), 1);
    });
}

BOOST_AUTO_TEST_CASE(test_clone_and_cast)
{
    BOOST_REQUIRE_NO_THROW({
        fwdpp::binomial_interval p(0, 1, 1e-3, false);
        auto c = p.clone();
        apply_callback(*c, rng.get());
    });
    fwdpp::binomial_interval p(0, 1, 1e-3, false);
    auto c = p.clone();
    BOOST_REQUIRE_EQUAL(c == nullptr, false);
    auto cast = dynamic_cast<fwdpp::binomial_interval*>(c.get());
    BOOST_REQUIRE_EQUAL(cast != nullptr, true);
    BOOST_REQUIRE_EQUAL(cast->beg, p.beg);
    BOOST_REQUIRE_EQUAL(cast->end, p.end);
    BOOST_REQUIRE_EQUAL(cast->prob, p.prob);
}

BOOST_AUTO_TEST_CASE(test_discrete_construction)
{
    BOOST_REQUIRE_NO_THROW(fwdpp::binomial_interval(0., 1., 1e-3, true));
    BOOST_REQUIRE_THROW(fwdpp::binomial_interval(0., 1., -1e-3, true),
                        std::invalid_argument);
    BOOST_REQUIRE_THROW(fwdpp::binomial_interval(1., 0., 1e-3, true),
                        std::invalid_argument);
    fwdpp::binomial_interval b(0., 1., 1e-3, true);
    BOOST_REQUIRE_EQUAL(b.discrete(), true);
}

BOOST_FIXTURE_TEST_CASE(test_discrete_call_operator, rng_fixture)
{
    fwdpp::binomial_interval p(0.0, 10.0, 1, true);
    auto rv = apply_callback(p, rng.get());
    BOOST_REQUIRE(!rv.empty());
    for (auto i : rv)
        {
            BOOST_REQUIRE_EQUAL(i, std::floor(i));
        }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(test_poisson_point, rng_fixture)

BOOST_AUTO_TEST_CASE(test_continuous_construction)
{
    BOOST_REQUIRE_NO_THROW(fwdpp::poisson_point(0., 1e-3, false));
    BOOST_REQUIRE_THROW(fwdpp::poisson_point(0., -1e-3, false), std::invalid_argument);
    BOOST_REQUIRE_THROW(
        fwdpp::poisson_point(1., std::numeric_limits<double>::quiet_NaN(), false),
        std::invalid_argument);
    fwdpp::poisson_point p(0., 1e-3, false);
    BOOST_REQUIRE_EQUAL(p.discrete(), false);
}

BOOST_AUTO_TEST_CASE(test_callback_nothrow)
{
    BOOST_REQUIRE_NO_THROW({
        fwdpp::poisson_point p(0, 1e-3, false);
        apply_callback(p, rng.get());
    });
}

BOOST_AUTO_TEST_CASE(test_clone_and_cast)
{
    BOOST_REQUIRE_NO_THROW({
        fwdpp::poisson_point p(0, 1e-3, false);
        auto c = p.clone();
        apply_callback(*c, rng.get());
    });
    fwdpp::poisson_point p(0, 1e-3, false);
    auto c = p.clone();
    BOOST_REQUIRE_EQUAL(c == nullptr, false);
    auto cast = dynamic_cast<fwdpp::poisson_point*>(c.get());
    BOOST_REQUIRE_EQUAL(cast != nullptr, true);
    BOOST_REQUIRE_EQUAL(cast->position, p.position);
    BOOST_REQUIRE_EQUAL(cast->mean, p.mean);
}

BOOST_AUTO_TEST_CASE(test_discrete_construction)
{
    BOOST_REQUIRE_NO_THROW(fwdpp::poisson_point(0., 1e-3, true));
    BOOST_REQUIRE_THROW(fwdpp::poisson_point(0., -1e-3, true), std::invalid_argument);
    BOOST_REQUIRE_THROW(
        fwdpp::poisson_point(1., std::numeric_limits<double>::quiet_NaN(), true),
        std::invalid_argument);
    fwdpp::poisson_point p(0., 1e-3, true);
    BOOST_REQUIRE_EQUAL(p.discrete(), true);
}

BOOST_FIXTURE_TEST_CASE(test_discrete_call_operator, rng_fixture)
{
    fwdpp::poisson_point p(0.6666, 10.0, true);
    auto rv = apply_callback(p, rng.get());
    BOOST_REQUIRE(!rv.empty());
    BOOST_REQUIRE_EQUAL(rv[0], 0.);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(test_binomial_point, rng_fixture)

BOOST_AUTO_TEST_CASE(test_continuous_construction)
{
    BOOST_REQUIRE_NO_THROW(fwdpp::binomial_point(0., 1e-3, false));
    BOOST_REQUIRE_THROW(fwdpp::binomial_point(0., -1e-3, false), std::invalid_argument);
    BOOST_REQUIRE_THROW(
        fwdpp::binomial_point(1., std::numeric_limits<double>::quiet_NaN(), false),
        std::invalid_argument);
    fwdpp::binomial_point b(0., 1e-3, false);
    BOOST_REQUIRE_EQUAL(b.discrete(), false);
}

BOOST_AUTO_TEST_CASE(test_callback_nothrow)
{
    BOOST_REQUIRE_NO_THROW({
        fwdpp::binomial_point p(0, 1, false);
        auto rv = apply_callback(p, rng.get());
        BOOST_REQUIRE_EQUAL(rv.size(), 1);
        BOOST_REQUIRE_EQUAL(rv.front(), 0);
    });
}

BOOST_AUTO_TEST_CASE(test_clone_and_cast)
{
    BOOST_REQUIRE_NO_THROW({
        fwdpp::binomial_point p(0, 1e-3, false);
        auto c = p.clone();
        apply_callback(*c, rng.get());
    });
    fwdpp::binomial_point p(0, 1e-3, false);
    auto c = p.clone();
    BOOST_REQUIRE_EQUAL(c == nullptr, false);
    auto cast = dynamic_cast<fwdpp::binomial_point*>(c.get());
    BOOST_REQUIRE_EQUAL(cast != nullptr, true);
    BOOST_REQUIRE_EQUAL(cast->position, p.position);
    BOOST_REQUIRE_EQUAL(cast->prob, p.prob);
}

BOOST_AUTO_TEST_CASE(test_discrete_construction)
{
    BOOST_REQUIRE_NO_THROW(fwdpp::binomial_point(0., 1e-3, true));
    BOOST_REQUIRE_THROW(fwdpp::binomial_point(0., -1e-3, true), std::invalid_argument);
    BOOST_REQUIRE_THROW(
        fwdpp::binomial_point(1., std::numeric_limits<double>::quiet_NaN(), true),
        std::invalid_argument);
    fwdpp::binomial_point b(0., 1e-3, true);
    BOOST_REQUIRE_EQUAL(b.discrete(), true);
}

BOOST_FIXTURE_TEST_CASE(test_discrete_call_operator, rng_fixture)
{
    fwdpp::binomial_point p(0.6666, 1.0, true);
    auto rv = apply_callback(p, rng.get());
    BOOST_REQUIRE(!rv.empty());
    BOOST_REQUIRE_EQUAL(rv[0], 0.);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(test_fixed_number_crossovers, rng_fixture)

BOOST_AUTO_TEST_CASE(test_continuous_construction)
{
    BOOST_REQUIRE_NO_THROW(fwdpp::fixed_number_crossovers(0, 1, 10, false));
    BOOST_REQUIRE_THROW(fwdpp::fixed_number_crossovers(1, 0, 10, false),
                        std::invalid_argument);
    BOOST_REQUIRE_THROW(fwdpp::fixed_number_crossovers(0, 1, -1, false),
                        std::invalid_argument);
    BOOST_REQUIRE_THROW(fwdpp::fixed_number_crossovers(
                            std::numeric_limits<double>::quiet_NaN(), 1, 1, false),
                        std::invalid_argument);
    BOOST_REQUIRE_THROW(fwdpp::fixed_number_crossovers(
                            1., std::numeric_limits<double>::quiet_NaN(), 1, false),
                        std::invalid_argument);
    fwdpp::fixed_number_crossovers f(0, 1, 10, false);
    BOOST_REQUIRE_EQUAL(f.discrete(), false);
}

BOOST_AUTO_TEST_CASE(test_callback_nothrow)
{
    BOOST_REQUIRE_NO_THROW({
        fwdpp::fixed_number_crossovers p(0, 1, 5, false);
        auto rv = apply_callback(p, rng.get());
        BOOST_REQUIRE_EQUAL(rv.size(), 5);
    });
}

BOOST_AUTO_TEST_CASE(test_clone_and_cast)
{
    BOOST_REQUIRE_NO_THROW({
        fwdpp::fixed_number_crossovers p(0, 1, 1, false);
        auto c = p.clone();
        apply_callback(*c, rng.get());
    });
    fwdpp::fixed_number_crossovers p(0, 1, 1, false);
    auto c = p.clone();
    BOOST_REQUIRE_EQUAL(c == nullptr, false);
    auto cast = dynamic_cast<fwdpp::fixed_number_crossovers*>(c.get());
    BOOST_REQUIRE_EQUAL(cast != nullptr, true);
    BOOST_REQUIRE_EQUAL(cast->beg, p.beg);
    BOOST_REQUIRE_EQUAL(cast->end, p.end);
    BOOST_REQUIRE_EQUAL(cast->nxovers, p.nxovers);
}

BOOST_AUTO_TEST_CASE(test_discrete_construction)
{
    BOOST_REQUIRE_NO_THROW(fwdpp::fixed_number_crossovers(0, 1, 10, true));
    BOOST_REQUIRE_THROW(fwdpp::fixed_number_crossovers(1, 0, 10, true),
                        std::invalid_argument);
    BOOST_REQUIRE_THROW(fwdpp::fixed_number_crossovers(0, 1, -1, true),
                        std::invalid_argument);
    BOOST_REQUIRE_THROW(fwdpp::fixed_number_crossovers(
                            std::numeric_limits<double>::quiet_NaN(), 1, 1, true),
                        std::invalid_argument);
    BOOST_REQUIRE_THROW(fwdpp::fixed_number_crossovers(
                            1., std::numeric_limits<double>::quiet_NaN(), 1, true),
                        std::invalid_argument);
    fwdpp::fixed_number_crossovers f(0, 1, 10, true);
    BOOST_REQUIRE_EQUAL(f.discrete(), true);
}

BOOST_FIXTURE_TEST_CASE(test_discrete_call_operator, rng_fixture)
{
    fwdpp::fixed_number_crossovers p(0, 1, 1, true);
    auto rv = apply_callback(p, rng.get());
    BOOST_REQUIRE_EQUAL(rv.size(), 1);
    BOOST_REQUIRE_EQUAL(rv[0], 0.);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(test_genetic_map, rng_fixture)

BOOST_AUTO_TEST_CASE(test_moving_vector)
{
    BOOST_REQUIRE_NO_THROW({
        std::vector<std::unique_ptr<fwdpp::genetic_map_unit>> callbacks;
        callbacks.emplace_back(new fwdpp::poisson_interval(0., 1., 1e-3, false));
        fwdpp::genetic_map gm(std::move(callbacks));
        BOOST_REQUIRE_EQUAL(gm.size(), 1);
        auto r = fwdpp::recbinder(std::cref(gm), rng.get());
        auto breakpoints = r();
    });
}

BOOST_AUTO_TEST_CASE(test_adding_poisson_interval)
{
    fwdpp::poisson_interval p(0, 1, 1e-1, false);
    fwdpp::genetic_map gm;
    gm.add_callback(p);
    BOOST_REQUIRE_EQUAL(gm.size(), 1);
    auto r = fwdpp::recbinder(std::cref(gm), rng.get());
    auto breakpoints = r();
}

BOOST_AUTO_TEST_CASE(test_adding_poisson_point)
{
    fwdpp::poisson_point p(0, 1, false);
    fwdpp::genetic_map gm;
    gm.add_callback(p);
    BOOST_REQUIRE_EQUAL(gm.size(), 1);
    auto r = fwdpp::recbinder(std::cref(gm), rng.get());
    auto breakpoints = r();
}

BOOST_AUTO_TEST_CASE(test_adding_binomial_point)
{
    fwdpp::binomial_point p(0, 1, false);
    fwdpp::genetic_map gm;
    gm.add_callback(p);
    BOOST_REQUIRE_EQUAL(gm.size(), 1);
    auto r = fwdpp::recbinder(std::cref(gm), rng.get());
    auto breakpoints = r();
}

BOOST_AUTO_TEST_SUITE_END()

