/*
  \file extensions.cc
  API checks on fwdpp's extensions sub-library.
*/

#include <config.h>
#include <boost/test/unit_test.hpp>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/extensions/regions.hpp>
#include <fwdpp/type_traits.hpp>
#include <limits>
#include "../fixtures/sugar_fixtures.hpp"

using namespace KTfwd;

using poptype = singlepop_popgenmut_fixture::poptype;
BOOST_FIXTURE_TEST_SUITE(test_extensions, singlepop_popgenmut_fixture)

// Check that extensions::discrete_mut_model::make_mut compiles
BOOST_AUTO_TEST_CASE(discrete_mut_model_test_1)
{
    // attempt
    extensions::discrete_mut_model dm({ 0, 1 }, { 1, 2 }, { 1, 0.5 }, {}, {},
                                      {}, {});
    //Check copy-constructible:
    decltype(dm) dm2(dm);
    //move-constructible:
    decltype(dm) dm3(std::move(dm2));
    auto rb = fwdpp_internal::make_mut_queue(pop.mcounts);
    KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);
    auto x = dm.make_mut(rb, pop.mutations, rng.get(), 0.001, 0., &generation,
                         pop.mut_lookup);
    static_assert(
        std::is_same<decltype(x), std::size_t>::value,
        "extensions::discrete_mut_model::make_muts must return a std::size_t");
}
// Check that extensions::discrete_mut_model::make_mut can be bound
BOOST_AUTO_TEST_CASE(discrete_mut_model_test_2)
{
    // attempt
    extensions::discrete_mut_model dm({ 0, 1 }, { 1, 2 }, { 1, 0.5 }, {}, {},
                                      {}, {});
    auto rb = fwdpp_internal::make_mut_queue(pop.mcounts);
    KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);
    auto mmodel = std::bind(
        &extensions::discrete_mut_model::
            make_mut<KTfwd::traits::recycling_bin_t<decltype(pop.mutations)>,
                     decltype(pop.mut_lookup), decltype(pop.mutations)>,
        &dm, rb, pop.mutations, rng.get(), 0.001, 0., &generation,
        std::ref(pop.mut_lookup));
    auto x = mmodel();
    static_assert(
        std::is_same<decltype(x), std::size_t>::value,
        "extensions::discrete_mut_model::make_muts must return a std::size_t");
}

// Check that extensions::discrete_mut_model::make_mut can be bound
// with placeholders, and that the resulting type is a valid
// mutation model
BOOST_AUTO_TEST_CASE(discrete_mut_model_test_3)
{
    // attempt
    extensions::discrete_mut_model dm({ 0, 1 }, { 1, 2 }, { 1, 0.5 }, {}, {},
                                      {}, {});
    auto rb = fwdpp_internal::make_mut_queue(pop.mcounts);
    KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);
    auto mmodel = std::bind(
        &extensions::discrete_mut_model::
            make_mut<KTfwd::traits::recycling_bin_t<decltype(pop.mutations)>,
                     decltype(pop.mut_lookup), decltype(pop.mutations)>,
        &dm, std::placeholders::_1, std::placeholders::_2, rng.get(), 0.001,
        0., &generation, std::ref(pop.mut_lookup));
    static_assert(
        traits::is_mutation_model<decltype(mmodel), poptype::mcont_t,
                                  poptype::gcont_t>::value,
        "error: type mutation_model is not a dispatchable mutation model "
        "type!");
    auto x = mmodel(rb, pop.mutations);
    static_assert(
        std::is_same<decltype(x), std::size_t>::value,
        "extensions::discrete_mut_model::make_muts must return a std::size_t");
}

// check return type of extensions::discrete_rec_model
BOOST_AUTO_TEST_CASE(discrete_rec_model_test_1)
{
    extensions::discrete_rec_model drm({ 0, 1 }, { 1, 2 }, { 1, 2 });
    KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);
    // use really big recombination rate here to ensure that return value is
    // not empty
    auto x
        = drm(rng.get(), 10.0, pop.gametes[0], pop.gametes[0], pop.mutations);
    static_assert(std::is_same<decltype(x), std::vector<double>>::value,
                  "extensions::dicrete_rec_model::operator() must return "
                  "std::vector<double>");
    BOOST_REQUIRE(x.empty()
                  || (x.back() == std::numeric_limits<double>::max()));
}

// test binding of extensions::discrete_rec_model::operator()
BOOST_AUTO_TEST_CASE(discrete_rec_model_test_2)
{
    extensions::discrete_rec_model drm({ 0, 1 }, { 1, 2 }, { 1, 2 });
    KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);
    auto bound = std::bind(
        &extensions::discrete_rec_model::
        operator()<poptype::gamete_t, decltype(pop.mutations)>,
        &drm, rng.get(), 0.001, pop.gametes[0], pop.gametes[0], pop.mutations);
    auto x = bound();
    static_assert(std::is_same<decltype(x), std::vector<double>>::value,
                  "extensions::dicrete_rec_model::operator() must return "
                  "std::vector<double>");
}

// test binding of extensions::discrete_rec_model::operator()
BOOST_AUTO_TEST_CASE(discrete_rec_model_test_3)
{
    extensions::discrete_rec_model drm({ 0, 1 }, { 1, 2 }, { 1, 2 });
    KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);
    auto bound
        = std::bind(&extensions::discrete_rec_model::
                    operator()<poptype::gamete_t, decltype(pop.mutations)>,
                    &drm, rng.get(), 0.001, std::placeholders::_1,
                    std::placeholders::_2, std::placeholders::_3);
    static_assert(traits::is_rec_model<decltype(bound), poptype::gamete_t,
                                       poptype::mcont_t>::value,
                  "extensions::discrete_rec_model::operator() is not a valid "
                  "recombination policy");
    auto x = bound(pop.gametes[0], pop.gametes[0], pop.mutations);
    static_assert(std::is_same<decltype(x), std::vector<double>>::value,
                  "extensions::dicrete_rec_model::operator() must return "
                  "std::vector<double>");
}

BOOST_AUTO_TEST_CASE(bound_drm_is_recmodel)
{
    extensions::discrete_rec_model drm({ 0, 1 }, { 1, 2 }, { 1, 1 });
    auto bound = extensions::bind_drm(drm, pop.gametes, pop.mutations,
                                      rng.get(), 0.001);
    static_assert(
        KTfwd::traits::
            is_rec_model<decltype(bound),
                         singlepop_popgenmut_fixture::poptype::gamete_t,
                         singlepop_popgenmut_fixture::poptype::mcont_t>::value,
        "bound object must be valid recombination model");
}

// Put it all together into a call to KTfwd::sample_diploid
/*
BOOST_AUTO_TEST_CASE( discrete_rec_model_test_4 )
{
  KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);

  extensions::discrete_mut_model dm({0,1},
                                    {1,2},
                                    {1,0.5},
                                    {},
                                    {},
                                    {},
                                    {}
                                    );

  auto mmodel =
std::bind(&extensions::discrete_mut_model::make_mut<KTfwd::traits::recycling_bin_t<decltype(pop.mutations)>,decltype(pop.mut_lookup),decltype(pop.mutations)>,
                          &dm,std::placeholders::_1,std::placeholders::_2,rng.get(),0.001,0.,0u,std::ref(pop.mut_lookup));

  extensions::discrete_rec_model drm( {0,1},
                                      {1,2},
                                      {1,2}
                                      );

  auto bound =
std::bind(&extensions::discrete_rec_model::operator()<poptype::gamete_t,decltype(pop.mutations)>,
                         &drm,
                         rng.get(),0.001,std::placeholders::_1,std::placeholders::_2,std::placeholders::_3);
  static_assert(
traits::is_rec_model<decltype(bound),poptype::gamete_t,poptype::mcont_t>::value,
                 "extensions::discrete_rec_model::operator() is not a valid
recombination policy" );
  auto x = bound(pop.gametes[0],pop.gametes[0],pop.mutations);
  static_assert(std::is_same<decltype(x),std::vector<double>>::value,
                "extensions::dicrete_rec_model::operator() must return
std::vector<double>");

  auto wbar = KTfwd::sample_diploid(rng.get(),
                                    pop.gametes,
                                    pop.diploids,
                                    pop.mutations,
                                    pop.mcounts,
                                    1000,
                                    0.001,
                                    mmodel,
                                    bound,
                                    std::bind(KTfwd::multiplicative_diploid(),std::placeholders::_1,std::placeholders::_2,
                                              std::placeholders::_3,2.),
                                    pop.neutral,
                                    pop.selected);
}
*/
// Put it all together into a call to KTfwd::sample_diploid,
// using both convenience fxns instead of the nasty templates
/*
BOOST_AUTO_TEST_CASE( discrete_rec_model_test_5 )
{
  KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);

  extensions::discrete_mut_model dm({0,1},
                                    {1,2},
                                    {1,0.5},
                                    {},
                                    {},
                                    {},
                                    {}
                                    );

  extensions::discrete_rec_model drm( {0,1},
                                      {1,2},
                                      {1,2}
                                      );

  auto wbar = KTfwd::sample_diploid(rng.get(),
                                    pop.gametes,
                                    pop.diploids,
                                    pop.mutations,
                                    pop.mcounts,
                                    1000,
                                    0.001,
                                    extensions::bind_dmm(dm,pop.mutations,pop.mut_lookup,
                                                         rng.get(),0.001,0.,0u),
                                    extensions::bind_drm(drm,pop.gametes,pop.mutations,
                                                         rng.get(),0.001),
                                    std::bind(KTfwd::multiplicative_diploid(),std::placeholders::_1,std::placeholders::_2,
                                              std::placeholders::_3,2.),
                                    pop.neutral,
                                    pop.selected);
}
*/
// Tests of raising exceptions
BOOST_AUTO_TEST_CASE(discrete_rec_model_constructor_should_throw)
{
    {
        BOOST_REQUIRE_THROW(
            extensions::discrete_rec_model drm({ 0 }, { 1, 2 }, { 1, 2 }),
            std::invalid_argument);
    }
    {
        BOOST_REQUIRE_THROW(
            extensions::discrete_rec_model drm({ 0, 1 }, { 1 }, { 1, 2 }),
            std::invalid_argument);
    }
    {
        BOOST_REQUIRE_THROW(extensions::discrete_rec_model drm(
                                { 0, 1 }, { 1, 2 }, { 1, 2, 3 }),
                            std::invalid_argument);
    }
}

BOOST_AUTO_TEST_CASE(discrete_mut_model_constructor_should_throw)
{
    {
        BOOST_REQUIRE_THROW(extensions::discrete_mut_model dm(
                                { 0, 1 }, { 1, 2 }, { 1 }, {}, {}, {}, {}),
                            std::invalid_argument);
    }
    {
        BOOST_REQUIRE_THROW(
            extensions::discrete_mut_model dm({ 0, 1 }, { 1, 2 }, { 1, 2 },
                                              // incorrect number of weights
                                              { 0, 1 }, { 1, 2 }, { 1 }, {}),
            std::invalid_argument);
    }
    {
        BOOST_REQUIRE_THROW(
            extensions::discrete_mut_model dm(
                { 0, 1 }, { 1, 2 }, { 1, 2 },
                // There are selected regions, but no "sh models"
                { 0, 1 }, { 1, 2 }, { 1 }, {}),
            std::invalid_argument);
    }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(unit_test_bind_vectors_regions,
                         multiloc_popgenmut_fixture)

BOOST_AUTO_TEST_CASE(bind_vec_drm_test_exceptions)
{
    {
        std::vector<double> recrates(3, 1e-4);

        BOOST_REQUIRE_THROW(
            auto bound = extensions::bind_vec_drm(
                vdrm, pop.gametes, pop.mutations, rng.get(), recrates),
            std::invalid_argument);
    }
}

BOOST_AUTO_TEST_SUITE_END()
