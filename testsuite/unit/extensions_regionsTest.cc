/*
  \file extensions.cc
  API checks on fwdpp's extensions sub-library.
*/

#include <config.h>
#include <boost/test/unit_test.hpp>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/extensions/regions.hpp>
#include <limits>
#include "../fixtures/sugar_fixtures.hpp"

using namespace KTfwd;

using poptype = singlepop<popgenmut>;
BOOST_FIXTURE_TEST_SUITE(test_extensions, singlepop_popgenmut_fixture)

// Check that extensions::discrete_mut_model::make_mut compiles
BOOST_AUTO_TEST_CASE(discrete_mut_model_test_1)
{
    // attempt
    extensions::discrete_mut_model dm({ 0, 1 }, { 1, 2 }, { 1, 0.5 }, {}, {},
                                      {}, {});
    auto rb = fwdpp_internal::make_mut_queue(pop.mcounts);
    KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);
    auto x = dm.make_mut(rb, pop.mutations, rng.get(), 0.001, 0., 0,
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
        &dm, rb, pop.mutations, rng.get(), 0.001, 0., 0u,
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
        0., 0u, std::ref(pop.mut_lookup));
    static_assert(
        traits::valid_mutation_model<decltype(mmodel), poptype::mcont_t,
                                     poptype::gcont_t>::value,
        "error: type mutation_model is not a dispatchable mutation model "
        "type!");
    auto x = mmodel(rb, pop.mutations);
    static_assert(
        std::is_same<decltype(x), std::size_t>::value,
        "extensions::discrete_mut_model::make_muts must return a std::size_t");
}

// Check that extensions::discrete_mut_model::make_mut can be bound
// with placeholders, that the resulting type is a valid
// mutation model, and can be passed to KTfwd::sample_diploid
/*
BOOST_AUTO_TEST_CASE( discrete_mut_model_test_4 )
{
  //attempt
  extensions::discrete_mut_model dm({0,1},
                                    {1,2},
                                    {1,0.5},
                                    {},
                                    {},
                                    {},
                                    {}
                                    );
  KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);
  auto mmodel =
std::bind(&extensions::discrete_mut_model::make_mut<KTfwd::traits::recycling_bin_t<decltype(pop.mutations)>,decltype(pop.mut_lookup),decltype(pop.mutations)>,
                          &dm,std::placeholders::_1,std::placeholders::_2,rng.get(),0.001,0.,0u,std::ref(pop.mut_lookup));
  static_assert(
traits::valid_mutation_model<decltype(mmodel),poptype::mcont_t,poptype::gcont_t>::value,
                 "error: type mutation_model is not a dispatchable mutation
model type!" );
  auto wbar = KTfwd::sample_diploid(rng.get(),
                                    pop.gametes,
                                    pop.diploids,
                                    pop.mutations,
                                    pop.mcounts,
                                    1000,
                                    0.001,
                                    mmodel,
                                    std::bind(KTfwd::poisson_xover(),rng.get(),0.001,0.,2.,
                                              std::placeholders::_1,std::placeholders::_2,std::placeholders::_3),
                                    std::bind(KTfwd::multiplicative_diploid(),std::placeholders::_1,std::placeholders::_2,
                                              std::placeholders::_3,2.),
                                    pop.neutral,
                                    pop.selected);
}
*/
/*
//Test the convenience fxn
BOOST_AUTO_TEST_CASE( discrete_mut_model_test_5 )
{
  //attempt
  extensions::discrete_mut_model dm({0,1},
                                    {1,2},
                                    {1,0.5},
                                    {},
                                    {},
                                    {},
                                    {}
                                    );
  KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);

  for(unsigned generation=0;generation<10000;++generation)
    {
      auto wbar = KTfwd::sample_diploid(rng.get(),
                                        pop.gametes,
                                        pop.diploids,
                                        pop.mutations,
                                        pop.mcounts,
                                        1000,
                                        0.001,
                                        extensions::bind_dmm(dm,pop.mutations,pop.mut_lookup,
                                                             rng.get(),0.001,0.,generation),
                                        std::bind(KTfwd::poisson_xover(),rng.get(),0.001,0.,2.,
                                                  std::placeholders::_1,std::placeholders::_2,std::placeholders::_3),
                                        std::bind(KTfwd::multiplicative_diploid(),std::placeholders::_1,std::placeholders::_2,
                                                  std::placeholders::_3,2.),
                                        pop.neutral,
                                        pop.selected);
      KTfwd::update_mutations(pop.mutations,pop.fixations,pop.fixation_times,pop.mut_lookup,pop.mcounts,generation,2000);
    }
  BOOST_REQUIRE_EQUAL(pop.mutations.size(),pop.mcounts.size());
  for(std::size_t i = 0 ; i < pop.mcounts.size() ; ++i )
    {
      if(pop.mcounts[i])
        {
          BOOST_REQUIRE( pop.mut_lookup.find(pop.mutations[i].pos) !=
pop.mut_lookup.end() );
        }
      else
        {
          BOOST_REQUIRE( pop.mut_lookup.find(pop.mutations[i].pos) ==
pop.mut_lookup.end() );
        }
    }

}
*/
/*
  Now, test discrete_mut_model's constructor that takes labels,
  a feature introduced in 0.4.9.  The purpose of this is to
  assign to mutation_base::xtra, which allows mutations to be integer-labelled.
*/
/*
BOOST_AUTO_TEST_CASE( discrete_mut_model_test_6 )
//This is an 'integration' test, I guess...
{
  //attempt
  extensions::discrete_mut_model dm({0,1},   //starts of 'neutral' regions
                                    {1,2},   //ends of 'neutral' regions
                                    {1,0.5}, //weights on 'neutral' regions
                                    {},      //starts of 'selected' regions
                                    {},      //stops of 'selected' regions
                                    {},      //weights on 'selected' regions
                                    {0,1},   //labels to put on mutations from
each of the 'neutral' regions
                                    {},      //labels to put on mutations from
each of the 'selected' regions
                                    {}       //vector of shmodels
                                    );

  //now, evolve the population
  KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);
  auto mmodel =
std::bind(&extensions::discrete_mut_model::make_mut<KTfwd::traits::recycling_bin_t<decltype(pop.mutations)>,decltype(pop.mut_lookup),decltype(pop.mutations)>,
                          &dm,std::placeholders::_1,std::placeholders::_2,rng.get(),0.001,0.,0u,std::ref(pop.mut_lookup));
  static_assert(
traits::valid_mutation_model<decltype(mmodel),poptype::mcont_t,poptype::gcont_t>::value,
                 "error: type mutation_model is not a dispatchable mutation
model type!" );
  auto wbar = KTfwd::sample_diploid(rng.get(),
                                    pop.gametes,
                                    pop.diploids,
                                    pop.mutations,
                                    pop.mcounts,
                                    1000,
                                    0.01,     //mutation rate--high so that
there are lots of mutations to test below...
                                    mmodel,
                                    std::bind(KTfwd::poisson_xover(),rng.get(),0.001,0.,2.,
                                              std::placeholders::_1,std::placeholders::_2,std::placeholders::_3),
                                    std::bind(KTfwd::multiplicative_diploid(),std::placeholders::_1,std::placeholders::_2,
                                              std::placeholders::_3,2.),
                                    pop.neutral,
                                    pop.selected);
  //Check that mutations in certain position intervals have the correct label
  for( const auto & m : pop.mutations )
    {
      if(m.pos < 1)
        {
          BOOST_REQUIRE_EQUAL(m.xtra,0);
        }
      else
        {
          BOOST_REQUIRE_EQUAL(m.xtra,1);
        }
    }
}
*/

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
    auto bound = std::bind(&extensions::discrete_rec_model::operator()
                               < poptype::gamete_t,
                           decltype(pop.mutations) >, &drm, rng.get(), 0.001,
                           pop.gametes[0], pop.gametes[0], pop.mutations);
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
    auto bound = std::bind(
        &extensions::discrete_rec_model::operator() < poptype::gamete_t,
        decltype(pop.mutations) >, &drm, rng.get(), 0.001,
        std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
    static_assert(traits::valid_rec_model<decltype(bound), poptype::gamete_t,
                                          poptype::mcont_t>::value,
                  "extensions::discrete_rec_model::operator() is not a valid "
                  "recombination policy");
    auto x = bound(pop.gametes[0], pop.gametes[0], pop.mutations);
    static_assert(std::is_same<decltype(x), std::vector<double>>::value,
                  "extensions::dicrete_rec_model::operator() must return "
                  "std::vector<double>");
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
traits::valid_rec_model<decltype(bound),poptype::gamete_t,poptype::mcont_t>::value,
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
            std::runtime_error);
    }
    {
        BOOST_REQUIRE_THROW(
            extensions::discrete_rec_model drm({ 0, 1 }, { 1 }, { 1, 2 }),
            std::runtime_error);
    }
    {
        BOOST_REQUIRE_THROW(extensions::discrete_rec_model drm(
                                { 0, 1 }, { 1, 2 }, { 1, 2, 3 }),
                            std::runtime_error);
    }
}

BOOST_AUTO_TEST_CASE(discrete_mut_model_constructor_should_throw)
{
    {
        BOOST_REQUIRE_THROW(extensions::discrete_mut_model dm(
                                { 0, 1 }, { 1, 2 }, { 1 }, {}, {}, {}, {}),
                            std::runtime_error);
    }
    {
        BOOST_REQUIRE_THROW(
            extensions::discrete_mut_model dm({ 0, 1 }, { 1, 2 }, { 1, 2 },
                                              // incorrect number of weights
                                              { 0, 1 }, { 1, 2 }, { 1 }, {}),
            std::runtime_error);
    }
    {
        BOOST_REQUIRE_THROW(
            extensions::discrete_mut_model dm(
                { 0, 1 }, { 1, 2 }, { 1, 2 },
                // There are selected regions, but no "sh models"
                { 0, 1 }, { 1, 2 }, { 1 }, {}),
            std::runtime_error);
    }
}

BOOST_AUTO_TEST_SUITE_END()
