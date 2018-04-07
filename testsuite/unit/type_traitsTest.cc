/*!
  \file type_traitsTest.cc
  \ingroup unit
  \brief Testing fwdpp/type_traits.hpp

  These tests make sure that the type traits
  actually return what we expect them to.
*/

#include <config.h>
#include "../fixtures/fwdpp_fixtures.hpp"
#include <fwdpp/fitness_models.hpp>
#include <fwdpp/mutate_recombine.hpp>
#include <fwdpp/poisson_xover.hpp>
#include <fwdpp/recbinder.hpp>
#include <boost/test/unit_test.hpp>
#include <fwdpp/type_traits.hpp>
#include <gsl/gsl_rng.h>

struct trivial_custom_diploid_invalid
/*!
  \ingroup unit
  Fails to define typedefs first_type and second_type
  \note see trivial_custom_diploid_valid
*/
{
};

struct trivial_custom_diploid_valid
/*!
  \ingroup unit
*/
{
    using first_type = std::size_t;
    using second_type = std::size_t;
};

BOOST_FIXTURE_TEST_SUITE(test_type_traits, standard_empty_single_deme_fixture)

BOOST_AUTO_TEST_CASE(is_diploid_test)
{
    auto v = fwdpp::traits::is_diploid<
        std::pair<std::size_t, std::size_t>>::value;
    BOOST_REQUIRE_EQUAL(v, true);
    v = fwdpp::traits::is_custom_diploid<
        std::pair<std::size_t, std::size_t>>::value;
    BOOST_REQUIRE_EQUAL(v, false);
}

BOOST_AUTO_TEST_CASE(is_gamete_test)
{
    auto v = fwdpp::traits::is_gamete<fwdpp::gamete>::value;
    BOOST_REQUIRE_EQUAL(v, true);
    v = fwdpp::traits::is_gamete<gcont_t::value_type>::value;
    BOOST_REQUIRE_EQUAL(v, true);
    v = fwdpp::traits::is_gamete<mtype>::value;
    BOOST_REQUIRE_EQUAL(v, false);
}

BOOST_AUTO_TEST_CASE(is_custom_diploid_test)
{
    auto v = fwdpp::traits::is_custom_diploid<
        trivial_custom_diploid_invalid>::value;
    BOOST_REQUIRE_EQUAL(v, false);
    v = fwdpp::traits::is_custom_diploid<trivial_custom_diploid_valid>::value;
    BOOST_REQUIRE_EQUAL(v, true);
}

BOOST_AUTO_TEST_CASE(is_mmodel_test)
{
    auto mut_recycling_bin = fwdpp::fwdpp_internal::make_mut_queue(mcounts);
    std::size_t i = 0;
    std::vector<double> next_mut_pos{ 0.0, 0.1, 0.2 };
    auto mmodel = [&next_mut_pos, &i](decltype(mut_recycling_bin) &rbin,
                                      std::vector<mtype> &__mvector) {
        // mutations are all neutral
        return fwdpp::fwdpp_internal::recycle_mutation_helper(
            rbin, __mvector, next_mut_pos[i++], 0.);
    };
    auto v
        = std::is_convertible<decltype(mmodel),
                              fwdpp::traits::mutation_model<mcont_t>>::value;
    BOOST_REQUIRE_EQUAL(v, true);

    v = fwdpp::traits::is_mutation_model<decltype(mmodel), mcont_t,
                                         gcont_t>::value;
    BOOST_REQUIRE_EQUAL(v, true);
}

BOOST_AUTO_TEST_CASE(is_mmodel_gamete_test)
{
    // Fake a mutation model requiring a gamete.
    // The function itself is invalid in implementation,
    // but the purpose here is to test the signature.
    auto m = [](fwdpp::traits::recycling_bin_t<mcont_t> &,
                const gcont_t::value_type &,
                const mcont_t &) -> std::size_t { return 1; };
    auto v = std::is_convertible<
        decltype(m),
        fwdpp::traits::mutation_model_gamete<mcont_t, gcont_t>>::value;
    BOOST_REQUIRE_EQUAL(v, true);
}

// Mocks a custom diploid
struct fake_dip
{
    using first_type = std::size_t;
    using second_type = first_type;
    first_type first;
    second_type second;
};

BOOST_AUTO_TEST_CASE(is_mmodel_diploid_test)
{
    // Fake a mutation model requiring a diploid.
    // The function itself is invalid in implementation,
    // but the purpose here is to test the signature.
    auto m = [](fwdpp::traits::recycling_bin_t<mcont_t> &, const fake_dip &,
                const gcont_t::value_type &,
                const mcont_t &) -> std::size_t { return 1; };
    auto v = std::is_convertible<
        decltype(m), fwdpp::traits::mutation_model_diploid<fake_dip, mcont_t,
                                                           gcont_t>>::value;
    BOOST_REQUIRE_EQUAL(v, true);
}

BOOST_AUTO_TEST_CASE(is_standard_fitness_model_test)
{
    auto fp = fwdpp::multiplicative_diploid(2.);
    auto v = std::is_convertible<
        decltype(fp),
        fwdpp::traits::fitness_fxn_t<dipvector_t, gcont_t, mcont_t>>::value;
    BOOST_REQUIRE_EQUAL(v, true);
}

BOOST_AUTO_TEST_CASE(is_not_fitness_model)
// These tests will simply fail to compile if they cannot pass.
// They are tests of compile-time concepts and not run-time
// expectations.
{
    auto v = fwdpp::traits::fitness_fxn<dipvector_t, std::vector<double>,
                                        mcont_t>();
    static_assert(std::is_void<decltype(v)::type>::value, "v must be void");
    auto ff = [](const dipvector_t &, const std::vector<double> &,
                 const mcont_t &) {};
    static_assert(
        !fwdpp::traits::is_fitness_fxn<decltype(ff), dipvector_t,
                                       std::vector<double>, mcont_t>::value,
        "foo");
}

BOOST_AUTO_TEST_CASE(is_empty_recmodel_test)
{
    const auto rm = fwdpp::recbinder(fwdpp::poisson_xover(1e-3, 0, 1), r);
    auto v = fwdpp::traits::is_rec_model<decltype(rm), dipvector_t::value_type,
                                         gcont_t::value_type, mcont_t>::value;
    BOOST_REQUIRE_EQUAL(v, true);
    // Now, we test that the types can be dispatched
    auto r1 = fwdpp::dispatch_recombination_policy(
        rm, dipvector_t::value_type(), gcont_t::value_type(0),
        gcont_t::value_type(0), mcont_t());
    v = std::is_same<decltype(r1), std::vector<double>>::value;
    BOOST_REQUIRE_EQUAL(v, true);
}

BOOST_AUTO_TEST_CASE(is_gamete_recmodel_test)
{
    const auto rm = fwdpp::recbinder(fwdpp::poisson_xover(1e-3, 0, 1), r);
    auto mock_rec
        = [&rm](const gcont_t::value_type &, const gcont_t::value_type &,
                const mcont_t &) -> decltype(rm()) { return rm(); };

    auto v = fwdpp::traits::is_rec_model<decltype(mock_rec),
                                         dipvector_t::value_type,
                                         gcont_t::value_type, mcont_t>::value;
    BOOST_REQUIRE_EQUAL(v, true);
    auto r2 = fwdpp::dispatch_recombination_policy(
        mock_rec, dipvector_t::value_type(), gcont_t::value_type(0),
        gcont_t::value_type(0), mcont_t());
    v = std::is_same<decltype(r2), std::vector<double>>::value;
    BOOST_REQUIRE_EQUAL(v, true);
}

BOOST_AUTO_TEST_CASE(is_diploid_recmodel_test)
{
    const auto rm = fwdpp::recbinder(fwdpp::poisson_xover(1e-3, 0, 1), r);
    auto mock_rec_diploid
        = [&rm](const dipvector_t::value_type &, const gcont_t::value_type &,
                const gcont_t::value_type &,
                const mcont_t &) -> decltype(rm()) { return rm(); };

    auto v = fwdpp::traits::is_rec_model<decltype(mock_rec_diploid),
                                         dipvector_t::value_type,
                                         gcont_t::value_type, mcont_t>::value;
    BOOST_REQUIRE_EQUAL(v, true);
    auto r2 = fwdpp::dispatch_recombination_policy(
        mock_rec_diploid, dipvector_t::value_type(), gcont_t::value_type(0),
        gcont_t::value_type(0), mcont_t());
    v = std::is_same<decltype(r2), std::vector<double>>::value;
}

BOOST_AUTO_TEST_CASE(is_not_recmodel_test)
{
    const auto rm = fwdpp::recbinder(fwdpp::poisson_xover(1e-3, 0, 1), r);
    // Test that this is not a valid rec policy
    auto mock_not_rec
        = [&rm](const int, int, const mcont_t &) -> decltype(rm()) {
        return rm();
    };
    auto v = fwdpp::traits::is_rec_model<
        decltype(mock_not_rec), dipvector_t::value_type, int, mcont_t>::value;
    BOOST_REQUIRE_EQUAL(v, false);
}

BOOST_AUTO_TEST_SUITE_END()
