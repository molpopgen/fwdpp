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

//Check that extensions::discrete_mut_model::make_mut compiles
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
//Check that extensions::discrete_mut_model::make_mut can be bound
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
  auto mmodel =  std::bind(&extensions::discrete_mut_model::make_mut<KTfwd::traits::recycling_bin_t<decltype(pop.mutations)>::type,decltype(pop.mut_lookup),decltype(pop.mutations)>,
			   &dm,rng.get(),0.001,0.,0u,rb,pop.mutations,std::ref(pop.mut_lookup));
  auto x = mmodel();
  static_assert( std::is_same<decltype(x),std::size_t>::value,
		 "extensions::discrete_mut_model::make_muts must return a std::size_t" );
}

//Check that extensions::discrete_mut_model::make_mut can be bound
//with placeholders, and that the resulting type is a valid
//mutation model
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
  auto mmodel = std::bind(&extensions::discrete_mut_model::make_mut<KTfwd::traits::recycling_bin_t<decltype(pop.mutations)>::type,decltype(pop.mut_lookup),decltype(pop.mutations)>,
			  &dm,rng.get(),0.001,0.,0u,std::placeholders::_1,std::placeholders::_2,std::ref(pop.mut_lookup));
  static_assert( traits::valid_mutation_model<decltype(mmodel),poptype::mcont_t,poptype::gcont_t>::value,
		 "error: type mutation_model is not a dispatchable mutation model type!" );
  auto x = mmodel(rb,pop.mutations);
  static_assert( std::is_same<decltype(x),std::size_t>::value,
		 "extensions::discrete_mut_model::make_muts must return a std::size_t" );
}

//Check that extensions::discrete_mut_model::make_mut can be bound
//with placeholders, that the resulting type is a valid
//mutation model, and can be passed to KTfwd::sample_diploid
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
  auto mmodel = std::bind(&extensions::discrete_mut_model::make_mut<KTfwd::traits::recycling_bin_t<decltype(pop.mutations)>::type,decltype(pop.mut_lookup),decltype(pop.mutations)>,
			  &dm,rng.get(),0.001,0.,0u,std::placeholders::_1,std::placeholders::_2,std::ref(pop.mut_lookup));
  static_assert( traits::valid_mutation_model<decltype(mmodel),poptype::mcont_t,poptype::gcont_t>::value,
		 "error: type mutation_model is not a dispatchable mutation model type!" );
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

//check return type of extensions::discrete_rec_model
BOOST_AUTO_TEST_CASE( discrete_rec_model_test_1 )
{
  poptype pop(1000);

  extensions::discrete_rec_model drm( {0,1},
				      {1,2},
				      {1,2}
				      );
  KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);
  auto x = drm(rng.get(),0.001,pop.gametes[0],pop.gametes[0],pop.mutations);
  static_assert(std::is_same<decltype(x),std::vector<double>>::value,
		"extensions::dicrete_rec_model::operator() must return std::vector<double>");
}

//test bindind of extensions::discrete_rec_model::operator()
BOOST_AUTO_TEST_CASE( discrete_rec_model_test_2 )
{
  poptype pop(1000);

  extensions::discrete_rec_model drm( {0,1},
				      {1,2},
				      {1,2}
				      );
  KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);
  auto bound = std::bind(&extensions::discrete_rec_model::operator()<poptype::gamete_t,decltype(pop.mutations)>,
			 &drm,
			 rng.get(),0.001,pop.gametes[0],pop.gametes[0],pop.mutations);	     
  auto x = bound();
  static_assert(std::is_same<decltype(x),std::vector<double>>::value,
		"extensions::dicrete_rec_model::operator() must return std::vector<double>");
}

//test binding of extensions::discrete_rec_model::operator()
BOOST_AUTO_TEST_CASE( discrete_rec_model_test_3 )
{
  poptype pop(1000);

  extensions::discrete_rec_model drm( {0,1},
				      {1,2},
				      {1,2}
				      );
  KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);
  auto bound = std::bind(&extensions::discrete_rec_model::operator()<poptype::gamete_t,decltype(pop.mutations)>,
			 &drm,
			 rng.get(),0.001,std::placeholders::_1,std::placeholders::_2,std::placeholders::_3);
  static_assert( traits::valid_rec_model<decltype(bound),poptype::gamete_t,poptype::mcont_t>::value,
		 "extensions::discrete_rec_model::operator() is not a valid recombination policy" );
  auto x = bound(pop.gametes[0],pop.gametes[0],pop.mutations);
  static_assert(std::is_same<decltype(x),std::vector<double>>::value,
		"extensions::dicrete_rec_model::operator() must return std::vector<double>");
}

//test binding of extensions::discrete_rec_model::operator()
BOOST_AUTO_TEST_CASE( discrete_rec_model_test_4 )
{
  poptype pop(1000);

  KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);
  
  extensions::discrete_mut_model dm({0,1},
				    {1,2},
				    {1,0.5},
				    {},
				    {},
				    {},
				    {}
				    );
  
  auto mmodel = std::bind(&extensions::discrete_mut_model::make_mut<KTfwd::traits::recycling_bin_t<decltype(pop.mutations)>::type,decltype(pop.mut_lookup),decltype(pop.mutations)>,
			  &dm,rng.get(),0.001,0.,0u,std::placeholders::_1,std::placeholders::_2,std::ref(pop.mut_lookup));
  
  extensions::discrete_rec_model drm( {0,1},
				      {1,2},
				      {1,2}
				      );

  auto bound = std::bind(&extensions::discrete_rec_model::operator()<poptype::gamete_t,decltype(pop.mutations)>,
			 &drm,
			 rng.get(),0.001,std::placeholders::_1,std::placeholders::_2,std::placeholders::_3);
  static_assert( traits::valid_rec_model<decltype(bound),poptype::gamete_t,poptype::mcont_t>::value,
		 "extensions::discrete_rec_model::operator() is not a valid recombination policy" );
  auto x = bound(pop.gametes[0],pop.gametes[0],pop.mutations);
  static_assert(std::is_same<decltype(x),std::vector<double>>::value,
		"extensions::dicrete_rec_model::operator() must return std::vector<double>");

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
