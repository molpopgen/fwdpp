/*!
  \file test_general_rec_variation.cc
  \ingroup unit
  \brief Unit tests for fwdpp::general_rec_variation.  Also tests the built-in
  policies accompanying this type.
*/

#include <boost/test/unit_test.hpp>
#include <fwdpp/general_rec_variation.hpp>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/mutate_recombine.hpp>
#include <config.h>

struct fixture
{
    const gsl_rng* r;
    const double not_a_num, inf;
    fwdpp::general_rec_variation recvar;
    fwdpp::poisson_interval pi;
    fwdpp::crossover_point ci1;
    fwdpp::crossover_point ci2;
    fwdpp::crossover_point ci3;
    fwdpp::haploid_genome g;
    std::pair<std::size_t, std::size_t> diploid;
    std::vector<fwdpp::mutation> mutations;

    fixture()
        : r{ gsl_rng_alloc(gsl_rng_taus) },
          not_a_num{ std::numeric_limits<double>::quiet_NaN() },
          inf{ std::numeric_limits<double>::infinity() }, recvar{},
          pi{ 1e-3, 0., 1. }, ci1{ 1e-3, 0.5 }, ci2{ 1e-3, 0.5, true },
          ci3{ 1e-3, 0.5, false }, g{ 0 }, diploid{ 0, 0 }, mutations{}
    {
        gsl_rng_set(r, 42);
    }
};

BOOST_FIXTURE_TEST_SUITE(test_general_rec_variation, fixture)

BOOST_AUTO_TEST_CASE(test_poisson_interval)
{
    BOOST_REQUIRE_THROW(fwdpp::poisson_interval(not_a_num, 0., 1.),
                        std::invalid_argument);
    BOOST_REQUIRE_THROW(fwdpp::poisson_interval(1e-3, not_a_num, 1.),
                        std::invalid_argument);
    BOOST_REQUIRE_THROW(fwdpp::poisson_interval(1e-3, 0., not_a_num),
                        std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(test_crossover_point)
{
    BOOST_REQUIRE_THROW(fwdpp::crossover_point(inf, -2),
                        std::invalid_argument);
    BOOST_REQUIRE_THROW(fwdpp::crossover_point(-2, not_a_num),
                        std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(test_basic_use)
{
    recvar.recmap.push_back(pi);
    recvar.recmap.push_back(ci1);
    recvar.recmap.push_back(ci2);
    recvar.recmap.push_back(ci3);
    auto x = recvar(r);
}

BOOST_AUTO_TEST_CASE(test_dispatch)
{
    const auto bound = [this]() { return recvar(r); };
    auto x = fwdpp::dispatch_recombination_policy(bound, std::cref(diploid),
                                                  std::cref(g), std::cref(g),
                                                  std::cref(mutations));
}

BOOST_AUTO_TEST_SUITE_END()
