/*
  \file extensions.cc
  API checks on fwdpp's extensions sub-library.
*/

#define BOOST_TEST_MODULE extensions_test
#define BOOST_TEST_DYN_LINK 

#include <config.h>
#include <iostream>
#include <cmath>
#include <limits>
#include <type_traits>
#include <boost/test/unit_test.hpp>
#include <fwdpp/diploid.hh>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/sugar/popgenmut.hpp>
#include <fwdpp/sugar/singlepop.hpp>
#include <fwdpp/extensions/regions.hpp>
#include <fwdpp/extensions/callbacks.hpp>
#include <limits>

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
  auto x = dm.make_mut(rb,pop.mutations,rng.get(),0.001,0.,
		       0,pop.mut_lookup);
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
  auto mmodel =  std::bind(&extensions::discrete_mut_model::make_mut<KTfwd::traits::recycling_bin_t<decltype(pop.mutations)>,decltype(pop.mut_lookup),decltype(pop.mutations)>,
			   &dm,rb,pop.mutations,rng.get(),0.001,0.,0u,std::ref(pop.mut_lookup));
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
  auto mmodel = std::bind(&extensions::discrete_mut_model::make_mut<KTfwd::traits::recycling_bin_t<decltype(pop.mutations)>,decltype(pop.mut_lookup),decltype(pop.mutations)>,
			  &dm,std::placeholders::_1,std::placeholders::_2,rng.get(),0.001,0.,0u,std::ref(pop.mut_lookup));
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
  auto mmodel = std::bind(&extensions::discrete_mut_model::make_mut<KTfwd::traits::recycling_bin_t<decltype(pop.mutations)>,decltype(pop.mut_lookup),decltype(pop.mutations)>,
			  &dm,std::placeholders::_1,std::placeholders::_2,rng.get(),0.001,0.,0u,std::ref(pop.mut_lookup));
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

//Test the convenience fxn
BOOST_AUTO_TEST_CASE( discrete_mut_model_test_5 )
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
	  BOOST_REQUIRE( pop.mut_lookup.find(pop.mutations[i].pos) != pop.mut_lookup.end() );
	}
      else
	{
	  BOOST_REQUIRE( pop.mut_lookup.find(pop.mutations[i].pos) == pop.mut_lookup.end() );
	}
    }
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

//Put it all together into a call to KTfwd::sample_diploid
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
  
  auto mmodel = std::bind(&extensions::discrete_mut_model::make_mut<KTfwd::traits::recycling_bin_t<decltype(pop.mutations)>,decltype(pop.mut_lookup),decltype(pop.mutations)>,
			  &dm,std::placeholders::_1,std::placeholders::_2,rng.get(),0.001,0.,0u,std::ref(pop.mut_lookup));
  
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

//Put it all together into a call to KTfwd::sample_diploid,
//using both convenience fxns instead of the nasty templates
BOOST_AUTO_TEST_CASE( discrete_rec_model_test_5 )
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
//Tests of fwdpp/extensions/callbacks.hpp

//This test makes sure that each type of callback compiles
BOOST_AUTO_TEST_CASE( vector_shmodel )
{
  //PS, uniform initialization rocks...
  std::vector< extensions::shmodel > callbacks {
    {extensions::constant(1.),extensions::constant(0.)},
      {extensions::exponential(1.),extensions::exponential(1.)},
	{extensions::uniform(1.,2.),extensions::uniform(1.,2.)},
	  {extensions::beta(1.,2.),extensions::beta(1.,2.)},             //defaults to factor = 1
	    {extensions::beta(1.,2.,0.25),extensions::beta(1.,2.,0.25)}, //pass all 3 params to constructor
	      {extensions::gaussian(1.),extensions::gaussian(1.)},
		{extensions::gamma(1.,0.1),extensions::gamma(1.,0.1)}
  };
}

//The callbacks can throw exceptions if their parameters aren't valid

BOOST_AUTO_TEST_CASE( callback_exceptions )
{
  {
    //inf
    BOOST_REQUIRE_THROW( extensions::constant(1./0.), 
			 std::runtime_error
			 );
  }
  
  {
    //nan
    BOOST_REQUIRE_THROW( extensions::constant(std::nan("")), 
     			 std::runtime_error
     			 );
  }

  {
    //first arg not finite
    BOOST_REQUIRE_THROW( extensions::uniform(1./0., 1.),
			 std::runtime_error
			 );
  }

  {
    //2nd arg not finite
    BOOST_REQUIRE_THROW( extensions::uniform(1., 1./0.),
			 std::runtime_error
			 );
  }

  {
    //min > max
    BOOST_REQUIRE_THROW( extensions::uniform(1., 0.99),
			 std::runtime_error
			 );
  }

  {
    //a not finite
    BOOST_REQUIRE_THROW( extensions::beta(std::nan(""),1.),
			 std::runtime_error
			 );
  }

  {
    //b not finite
    BOOST_REQUIRE_THROW( extensions::beta(1.,std::nan("")),
			 std::runtime_error
			 );
  }

  {
    //f not finite
    BOOST_REQUIRE_THROW( extensions::beta(1.,1.,std::nan("")),
			 std::runtime_error
			 );
  }

  {
    //a <= 0.
    BOOST_REQUIRE_THROW( extensions::beta(0.,1.),
			 std::runtime_error
			 );
  }

  {
    //b <= 0.
    BOOST_REQUIRE_THROW( extensions::beta(1.,0.),
			 std::runtime_error
			 );
  }

  {
    //f <= 0.
    BOOST_REQUIRE_THROW( extensions::beta(1.,1.,0.),
			 std::runtime_error
			 );
  }

  {
    //sd = 0
    BOOST_REQUIRE_THROW( extensions::gaussian(0.),
			 std::runtime_error
			 );
  }
  
  {
    //sd < 0
    BOOST_REQUIRE_THROW( extensions::gaussian(-1e-6),
			 std::runtime_error
			 );
  }

  {
    //sd not finite
    BOOST_REQUIRE_THROW( extensions::gaussian(1./0.),
			 std::runtime_error
			 );
  }

  {
    //mean not finite
    BOOST_REQUIRE_THROW( extensions::gamma(1./0.,1.),
			 std::runtime_error
			 );
  }

  {
    //shape not finite
    BOOST_REQUIRE_THROW( extensions::gamma(1.0,1./0.),
			 std::runtime_error
			 );
  }

  {
    //!(shape>0)
    BOOST_REQUIRE_THROW( extensions::gamma(1.0,0.),
			 std::runtime_error
			 );
  }
}
