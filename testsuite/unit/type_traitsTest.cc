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
#include <fwdpp/recombination.hpp>
#include <boost/test/unit_test.hpp>
#include <fwdpp/type_traits.hpp>
#include <gsl/gsl_rng.h>

struct trivial_custom_diploid_invalid : public KTfwd::tags::custom_diploid_t
/*!
  \ingroup unit
  Fails to define typedefs first_type and second_type
  \note see trivial_custom_diploid_valid
*/
{
};

struct trivial_custom_diploid_valid : public KTfwd::tags::custom_diploid_t
/*!
  \ingroup unit
*/
{
    using first_type = std::size_t;
    using second_type = std::size_t;
};

BOOST_FIXTURE_TEST_SUITE(test_type_traits, standard_empty_single_deme_fixture)

BOOST_AUTO_TEST_CASE(is_diploid_like_test)
{
    auto v = KTfwd::traits::is_diploid_like<std::pair<std::size_t,
                                                      std::size_t>>::value;
    BOOST_REQUIRE_EQUAL(v, true);
    v = KTfwd::traits::is_custom_diploid_t<std::pair<std::size_t,
                                                     std::size_t>>::value;
    BOOST_REQUIRE_EQUAL(v, false);
}

BOOST_AUTO_TEST_CASE(is_gamete_test)
{
    auto v = KTfwd::traits::is_gamete_t<KTfwd::gamete>::value;
    BOOST_REQUIRE_EQUAL(v, true);
    v = KTfwd::traits::is_gamete_t<gcont_t::value_type>::value;
    BOOST_REQUIRE_EQUAL(v, true);
    v = KTfwd::traits::is_gamete_t<mtype>::value;
    BOOST_REQUIRE_EQUAL(v, false);
}

BOOST_AUTO_TEST_CASE(is_custom_diploid_test)
{
    auto v = KTfwd::traits::
        is_custom_diploid_t<trivial_custom_diploid_invalid>::value;
    BOOST_REQUIRE_EQUAL(v, false);
    v = KTfwd::traits::is_custom_diploid_t<trivial_custom_diploid_valid>::
        value;
    BOOST_REQUIRE_EQUAL(v, true);
}

BOOST_AUTO_TEST_CASE(is_mmodel_test)
{
    auto mut_recycling_bin = KTfwd::fwdpp_internal::make_mut_queue(mcounts);
    std::size_t i = 0;
    std::vector<double> next_mut_pos{ 0.0, 0.1, 0.2 };
    auto mmodel = [&next_mut_pos, &i](decltype(mut_recycling_bin) &rbin,
                                      std::vector<mtype> &__mvector) {
        // mutations are all neutral
        return KTfwd::fwdpp_internal::recycle_mutation_helper(
            rbin, __mvector, next_mut_pos[i++], 0.);
    };
    auto v = std::is_convertible<decltype(mmodel),
                                 KTfwd::traits::mmodel_t<mcont_t>>::value;
    BOOST_REQUIRE_EQUAL(v, true);

    v = KTfwd::traits::valid_mutation_model<decltype(mmodel), mcont_t,
                                            gcont_t>::value;
    BOOST_REQUIRE_EQUAL(v, true);
}

BOOST_AUTO_TEST_CASE(is_standard_fitness_model_test)
{
    auto fp = std::bind(KTfwd::multiplicative_diploid(), std::placeholders::_1,
                        std::placeholders::_2, std::placeholders::_3, 2.);
    auto v = std::is_convertible<decltype(fp),
                                 KTfwd::traits::fitness_fxn_t<dipvector_t,
                                                              gcont_t,
                                                              mcont_t>>::value;
    BOOST_REQUIRE_EQUAL(v, true);
}

BOOST_AUTO_TEST_CASE(is_recmodel_test)
{
    auto rm = std::bind(KTfwd::poisson_xover(), r, 1e-2, 0., 1.,
                        std::placeholders::_1, std::placeholders::_2,
                        std::placeholders::_3);
    // auto v =
    // std::is_convertible<decltype(rm),KTfwd::traits::recmodel_t<singlepop_t::gcont_t,singlepop_t::mcont_t>
    // >::value;
    auto v = std::is_convertible<decltype(rm),
                                 KTfwd::traits::recmodel_t<gcont_t,
                                                           mcont_t>>::value;
    BOOST_REQUIRE_EQUAL(v, true);
    v = KTfwd::traits::valid_rec_model<decltype(rm), gcont_t, mcont_t>::value;
    BOOST_REQUIRE_EQUAL(v, true);
}

BOOST_AUTO_TEST_SUITE_END()
