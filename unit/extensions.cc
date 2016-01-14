/*
  \file extensions.cc
  API checks on fwdpp's extensions sub-library.

  This test passes if all blocks compile.
*/

#define BOOST_TEST_MODULE extensions_test
#define BOOST_TEST_DYN_LINK 

#include <config.h>
#include <iostream>
#include <type_traits>
#include <boost/test/unit_test.hpp>
#include <fwdpp/diploid.hh>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/sugar/popgenmut.hpp>
#include <fwdpp/sugar/singlepop.hpp>
#include <fwdpp/extensions/regions.hpp>

using namespace KTfwd;

using poptype = singlepop<popgenmut>;

BOOST_AUTO_TEST_CASE( discrete_mut_model_test_1 )
{
  poptype pop(1000);

  //attempt
  extensions::discrete_mut_model dm({0,1},
				    {1,2},
				    {1,0.5},
				    {},
				    {},
				    {},
				    {}
				    );
  auto rb = fwdpp_internal::make_mut_queue(pop.mcounts);
  KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);
  auto x = dm.make_mut(rng.get(),0.001,0.,
		       0,rb,pop.mutations,pop.mut_lookup);
  static_assert( std::is_same<decltype(x),std::size_t>::value,
		 "extensions::discrete_mut_model::make_muts must return a std::size_t" );
}

BOOST_AUTO_TEST_CASE( discrete_mut_model_test_2 )
{
  poptype pop(1000);

  //attempt
  extensions::discrete_mut_model dm({0,1},
				    {1,2},
				    {1,0.5},
				    {},
				    {},
				    {},
				    {}
				    );
  auto rb = fwdpp_internal::make_mut_queue(pop.mcounts);
  KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);
  auto mmodel =  std::bind(&extensions::discrete_mut_model::make_mut<KTfwd::recycling_bin_t<decltype(pop.mutations)>::type,decltype(pop.mut_lookup),decltype(pop.mutations)>,
			   &dm,rng.get(),0.001,0.,0u,rb,pop.mutations,std::ref(pop.mut_lookup));
  auto x = mmodel();
  static_assert( std::is_same<decltype(x),std::size_t>::value,
		 "extensions::discrete_mut_model::make_muts must return a std::size_t" );
}

BOOST_AUTO_TEST_CASE( discrete_mut_model_test_3 )
{
  poptype pop(1000);

  //attempt
  extensions::discrete_mut_model dm({0,1},
				    {1,2},
				    {1,0.5},
				    {},
				    {},
				    {},
				    {}
				    );
  auto rb = fwdpp_internal::make_mut_queue(pop.mcounts);
  KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);
  auto mmodel = std::bind(&extensions::discrete_mut_model::make_mut<KTfwd::recycling_bin_t<decltype(pop.mutations)>::type,decltype(pop.mut_lookup),decltype(pop.mutations)>,
			  &dm,rng.get(),0.001,0.,0u,std::placeholders::_1,std::placeholders::_2,std::ref(pop.mut_lookup));
  static_assert( traits::valid_mutation_model<decltype(mmodel),poptype::mcont_t,poptype::gcont_t>::value,
		 "error: type mutation_model is not a dispatchable mutation model type!" );
  auto x = mmodel(rb,pop.mutations);
  static_assert( std::is_same<decltype(x),std::size_t>::value,
		 "extensions::discrete_mut_model::make_muts must return a std::size_t" );
}

BOOST_AUTO_TEST_CASE( discrete_mut_model_test_4 )
{
  poptype pop(1000);

  //attempt
  extensions::discrete_mut_model dm({0,1},
				    {1,2},
				    {1,0.5},
				    {},
				    {},
				    {},
				    {}
				    );
  auto rb = fwdpp_internal::make_mut_queue(pop.mcounts);
  KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);
  auto mmodel = std::bind(&extensions::discrete_mut_model::make_mut<KTfwd::recycling_bin_t<decltype(pop.mutations)>::type,decltype(pop.mut_lookup),decltype(pop.mutations)>,
			  &dm,rng.get(),0.001,0.,0u,std::placeholders::_1,std::placeholders::_2,std::ref(pop.mut_lookup));
  static_assert( traits::valid_mutation_model<decltype(mmodel),poptype::mcont_t,poptype::gcont_t>::value,
		 "error: type mutation_model is not a dispatchable mutation model type!" );
  auto wbar = KTfwd::sample_diploid(rng.get(),
				    pop.gametes,  //non-const reference to gametes
				    pop.diploids, //non-const reference to diploids
				    pop.mutations, //non-const reference to mutations
				    pop.mcounts,
				    1000,     //current pop size, remains constant
				    0.001,    //mutation rate per gamete
				    mmodel,
				    std::bind(KTfwd::poisson_xover(),rng.get(),0.001,0.,2.,
					      std::placeholders::_1,std::placeholders::_2,std::placeholders::_3),
				    std::bind(KTfwd::multiplicative_diploid(),std::placeholders::_1,std::placeholders::_2,
					      std::placeholders::_3,2.),
				    pop.neutral,
				    pop.selected);
}
