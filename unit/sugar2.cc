//! \file sugar2.cc \ingroup unit
#define BOOST_TEST_MODULE sugarTest2
#define BOOST_TEST_DYN_LINK 

#include <config.h>
#include <boost/test/unit_test.hpp>
#include <gsl/gsl_rng.h>
#include <fwdpp/diploid.hh>
#define FWDPP_SUGAR_USE_BOOST
#include <fwdpp/sugar/singlepop.hpp>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/sugar/infsites.hpp>

using mutation_with_age = KTfwd::popgenmut;
using mwriter = KTfwd::mutation_writer;
using mreader = KTfwd::mutation_reader<mutation_with_age>;
using poptype = KTfwd::singlepop_serialized<mutation_with_age,mwriter,mreader>;

BOOST_AUTO_TEST_CASE( singlepop_serialized_copy_construct_test )
{
  poptype pop(1000);

  KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);
  
  //Evolve for 10 generations
  std::function<double(void)> recmap = std::bind(gsl_rng_uniform,rng);
  for( unsigned generation= 0 ; generation < 10 ; ++generation )
    {
      double wbar = KTfwd::sample_diploid(rng,
					  &pop.gametes,
					  &pop.diploids,
					  &pop.mutations,
					  1000,
					  0.005,
					  std::bind(KTfwd::infsites(),rng,std::placeholders::_1,&pop.mut_lookup,generation,
						    0.005,0.,[](gsl_rng * r){return gsl_rng_uniform(r);},[](gsl_rng * r){return 0.;},[](gsl_rng * r){return 0.;}),
					  std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,
						    &pop.gametes,
						    0., //no rec
						    rng,
						    recmap),
					  std::bind(KTfwd::insert_at_end<poptype::mutation_t,poptype::mlist_t>,std::placeholders::_1,std::placeholders::_2),
					  std::bind(KTfwd::insert_at_end<poptype::gamete_t,poptype::glist_t>,std::placeholders::_1,std::placeholders::_2),
					  std::bind(KTfwd::multiplicative_diploid(),std::placeholders::_1,std::placeholders::_2,2.),
					  std::bind(KTfwd::mutation_remover(),std::placeholders::_1,0,2*pop.N));
      KTfwd::remove_fixed_lost(&pop.mutations,&pop.fixations,&pop.fixation_times,&pop.mut_lookup,generation,2*pop.N);
    }

  //attempt to copy
  KTfwd::serialize s;
  poptype pop2(pop);

  BOOST_REQUIRE(pop.mutations.size() == pop2.mutations.size());
  BOOST_REQUIRE(pop.gametes.size() == pop2.gametes.size());
  BOOST_REQUIRE(pop.diploids.size() == pop2.diploids.size());
  BOOST_REQUIRE(pop.mut_lookup == pop2.mut_lookup);
  //Compare the mutations
  for( auto m1 = pop.mutations.begin(),m2 = pop2.mutations.begin() ; m1 != pop.mutations.end() ; ++m1,++m2 )
    {
      BOOST_CHECK_EQUAL( m1->pos, m2->pos );
      BOOST_CHECK_EQUAL( m1->n, m2->n );
    }
  
  //Compare the gametes
  for( auto g1 = pop.gametes.begin(),g2 = pop2.gametes.begin() ; g1 != pop.gametes.end() ; ++g1,++g2 )
    {
      BOOST_CHECK( g1 != g2 );
      for( auto m1 = g1->mutations.begin(),m2=g2->mutations.begin() ; m1 != g1->mutations.end() ; ++m1,++m2 )
	{
	  BOOST_CHECK( m1 != m2 );
	  BOOST_CHECK( std::distance(pop.mutations.begin(),*m1) == std::distance(pop2.mutations.begin(),*m2) );
	  BOOST_CHECK_EQUAL( (*m1)->pos, (*m2)->pos );
	  BOOST_CHECK_EQUAL( (*m1)->n, (*m2)->n ); 
	}
    }
  
  //Compare the diploids
  for( auto d1 = pop.diploids.begin(),d2 = pop2.diploids.begin() ; d1 != pop.diploids.end() ; ++d1,++d2 )
    {
      BOOST_CHECK(d1->first != d2->first);
      BOOST_CHECK(d1->second != d2->second);
      BOOST_CHECK( std::distance( pop.gametes.begin(),d1->first ) == std::distance( pop2.gametes.begin(),d2->first ) );
      BOOST_CHECK( std::distance( pop.gametes.begin(),d1->second ) == std::distance( pop2.gametes.begin(),d2->second ) );
      for( auto m1 = d1->first->mutations.begin(),m2=d2->first->mutations.begin() ; m1 != d1->first->mutations.end() ; ++m1,++m2 )
	{
	  BOOST_CHECK( m1 != m2 );
	  BOOST_CHECK( std::distance(pop.mutations.begin(),*m1) == std::distance(pop2.mutations.begin(),*m2) );
	  BOOST_CHECK_EQUAL( (*m1)->pos, (*m2)->pos );
	  BOOST_CHECK_EQUAL( (*m1)->n, (*m2)->n ); 
	}
      for( auto m1 = d1->second->mutations.begin(),m2=d2->second->mutations.begin() ; m1 != d1->second->mutations.end() ; ++m1,++m2 )
	{
	  BOOST_CHECK( m1 != m2 );
	  BOOST_CHECK( std::distance(pop.mutations.begin(),*m1) == std::distance(pop2.mutations.begin(),*m2) );
	  BOOST_CHECK_EQUAL( (*m1)->pos, (*m2)->pos );
	  BOOST_CHECK_EQUAL( (*m1)->n, (*m2)->n ); 
	}
    }
}

BOOST_AUTO_TEST_CASE( singlepop_sugar_assignment_test )
{
  poptype pop(1000);

  KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);
    
  //Evolve for 10 generations
  std::function<double(void)> recmap = std::bind(gsl_rng_uniform,rng);
  for( unsigned generation= 0 ; generation < 10 ; ++generation )
    {
      double wbar = KTfwd::sample_diploid(rng,
					  &pop.gametes,
					  &pop.diploids,
					  &pop.mutations,
					  1000,
					  0.005,
					  std::bind(KTfwd::infsites(),rng,std::placeholders::_1,&pop.mut_lookup,generation,
						    0.005,0.,[](gsl_rng * r){return gsl_rng_uniform(r);},[](gsl_rng * r){return 0.;},[](gsl_rng * r){return 0.;}),
					  std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,
						    &pop.gametes,
						    0., //no rec
						    rng,
						    recmap),
					  std::bind(KTfwd::insert_at_end<poptype::mutation_t,poptype::mlist_t>,std::placeholders::_1,std::placeholders::_2),
					  std::bind(KTfwd::insert_at_end<poptype::gamete_t,poptype::glist_t>,std::placeholders::_1,std::placeholders::_2),
					  std::bind(KTfwd::multiplicative_diploid(),std::placeholders::_1,std::placeholders::_2,2.),
					  std::bind(KTfwd::mutation_remover(),std::placeholders::_1,0,2*pop.N));
      KTfwd::remove_fixed_lost(&pop.mutations,&pop.fixations,&pop.fixation_times,&pop.mut_lookup,generation,2*pop.N);
    }

  //attempt to copy
  KTfwd::serialize s;
  poptype pop2 = pop;

  BOOST_REQUIRE(pop.mutations.size() == pop2.mutations.size());
  BOOST_REQUIRE(pop.gametes.size() == pop2.gametes.size());
  BOOST_REQUIRE(pop.diploids.size() == pop2.diploids.size());
  BOOST_REQUIRE(pop.mut_lookup == pop2.mut_lookup);
  //Compare the mutations
  for( auto m1 = pop.mutations.begin(),m2 = pop2.mutations.begin() ; m1 != pop.mutations.end() ; ++m1,++m2 )
    {
      BOOST_CHECK_EQUAL( m1->pos, m2->pos );
      BOOST_CHECK_EQUAL( m1->n, m2->n );
    }
  
  //Compare the gametes
  for( auto g1 = pop.gametes.begin(),g2 = pop2.gametes.begin() ; g1 != pop.gametes.end() ; ++g1,++g2 )
    {
      BOOST_CHECK( g1 != g2 );
      for( auto m1 = g1->mutations.begin(),m2=g2->mutations.begin() ; m1 != g1->mutations.end() ; ++m1,++m2 )
	{
	  BOOST_CHECK( m1 != m2 );
	  BOOST_CHECK( std::distance(pop.mutations.begin(),*m1) == std::distance(pop2.mutations.begin(),*m2) );
	  BOOST_CHECK_EQUAL( (*m1)->pos, (*m2)->pos );
	  BOOST_CHECK_EQUAL( (*m1)->n, (*m2)->n ); 
	}
    }
  
  //Compare the diploids
  for( auto d1 = pop.diploids.begin(),d2 = pop2.diploids.begin() ; d1 != pop.diploids.end() ; ++d1,++d2 )
    {
      BOOST_CHECK(d1->first != d2->first);
      BOOST_CHECK(d1->second != d2->second);
      BOOST_CHECK( std::distance( pop.gametes.begin(),d1->first ) == std::distance( pop2.gametes.begin(),d2->first ) );
      BOOST_CHECK( std::distance( pop.gametes.begin(),d1->second ) == std::distance( pop2.gametes.begin(),d2->second ) );
      for( auto m1 = d1->first->mutations.begin(),m2=d2->first->mutations.begin() ; m1 != d1->first->mutations.end() ; ++m1,++m2 )
	{
	  BOOST_CHECK( m1 != m2 );
	  BOOST_CHECK( std::distance(pop.mutations.begin(),*m1) == std::distance(pop2.mutations.begin(),*m2) );
	  BOOST_CHECK_EQUAL( (*m1)->pos, (*m2)->pos );
	  BOOST_CHECK_EQUAL( (*m1)->n, (*m2)->n ); 
	}
      for( auto m1 = d1->second->mutations.begin(),m2=d2->second->mutations.begin() ; m1 != d1->second->mutations.end() ; ++m1,++m2 )
	{
	  BOOST_CHECK( m1 != m2 );
	  BOOST_CHECK( std::distance(pop.mutations.begin(),*m1) == std::distance(pop2.mutations.begin(),*m2) );
	  BOOST_CHECK_EQUAL( (*m1)->pos, (*m2)->pos );
	  BOOST_CHECK_EQUAL( (*m1)->n, (*m2)->n ); 
	}
    }

  poptype pop3(std::move(pop));
  BOOST_REQUIRE(pop2.mut_lookup == pop3.mut_lookup);
  BOOST_REQUIRE(pop.mut_lookup != pop3.mut_lookup); //Assert that it moved
  BOOST_REQUIRE(pop.mut_lookup.empty());
}
