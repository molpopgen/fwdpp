/*!
  \file sugar3.cc
  \ingroup unit
  \brief Testing KTfwd::multiloc and KTfwd::multiloc_serialized
*/

#define BOOST_TEST_MODULE sugarTest3
#define BOOST_TEST_DYN_LINK 

#include <iostream>
#include <config.h>
#include <functional>
#include <algorithm>
#include <boost/test/unit_test.hpp>
#include <fwdpp/diploid.hh>
#ifndef USE_STANDARD_CONTAINERS //from config.h
#define FWDPP_SUGAR_USE_BOOST
#endif
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/sugar/multiloc.hpp>
#include <fwdpp/sugar/infsites.hpp>

using mwriter = KTfwd::mutation_writer;
using mreader = KTfwd::mutation_reader<KTfwd::popgenmut>;

//Fitness function
struct no_selection_multi
{
  typedef double result_type;
  template< typename diploid_type >
  inline double operator()( diploid_type & diploid ) const
  {
    return 1.;
  }
};

BOOST_AUTO_TEST_CASE( metapop_sugar_test1 )
{
  using poptype = KTfwd::multiloc<KTfwd::popgenmut>;
  KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);

  poptype pop(1000,4);
  BOOST_REQUIRE_EQUAL( pop.gametes[0].size(),1 );
  BOOST_REQUIRE_EQUAL( pop.gametes[1].size(),1 );
  BOOST_REQUIRE_EQUAL( pop.gametes[2].size(),1 );
  BOOST_REQUIRE_EQUAL( pop.gametes[3].size(),1 );
  BOOST_REQUIRE_EQUAL( pop.diploids.size(), 1000 );
  for( unsigned i = 0 ; i < pop.diploids.size() ; ++i )
    {
      BOOST_REQUIRE_EQUAL( pop.diploids[i].size(), 4 );
      auto gptr = pop.gametes.begin();
      for( unsigned j = 0 ; j < 4 ; ++j, ++gptr )
	{
	  BOOST_REQUIRE( pop.diploids[i][j].first == gptr->begin() );
	  BOOST_REQUIRE( pop.diploids[i][j].second == gptr->begin() );
	}
    }
  BOOST_REQUIRE( pop.mutations.empty() == true );
  unsigned generation = 0;
  
  //Mutation model for 4 loci
  std::vector< std::function<KTfwd::popgenmut(typename poptype::mlist_t *)> > mutmodels {
    //Locus 0: positions Uniform [0,1)
    std::bind(KTfwd::infsites(),rng,std::placeholders::_1,&pop.mut_lookup,&generation,
	      0.005,0.,[](gsl_rng * r){return gsl_rng_uniform(r);},[](gsl_rng * r){return 0.;},[](gsl_rng * r){return 0.;}) ,
      //Locus 1: positions Uniform [1,2)
      std::bind(KTfwd::infsites(),rng,std::placeholders::_1,&pop.mut_lookup,&generation,
      		0.005,0.,[](gsl_rng * r){return gsl_ran_flat(r,1.,2.);},[](gsl_rng * r){return 0.;},[](gsl_rng * r){return 0.;}),
      //Locus 2: positions Uniform [2,3)
      std::bind(KTfwd::infsites(),rng,std::placeholders::_1,&pop.mut_lookup,&generation,
      		0.005,0.,[](gsl_rng * r){return gsl_ran_flat(r,2.,3.);},[](gsl_rng * r){return 0.;},[](gsl_rng * r){return 0.;}),
      //Locus 3: positions Uniform [3,4)
      std::bind(KTfwd::infsites(),rng,std::placeholders::_1,&pop.mut_lookup,&generation,
      		0.005,0.,[](gsl_rng * r){return gsl_ran_flat(r,3.,4.);},[](gsl_rng * r){return 0.;},[](gsl_rng * r){return 0.;})
  };
  
  //Within-locus recombination models for 4 loci
  std::vector< std::function<unsigned(typename poptype::glist_t::iterator &,
				      typename poptype::glist_t::iterator &)> > recmodels {
    //Locus 0: positions = uniform [0,1)
    std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,&pop.gametes[0],0.005,rng,[&rng](){ return gsl_rng_uniform(rng); }),
      //Locus 1: positions = uniform [1,2)
      std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,&pop.gametes[1],0.005,rng,[&rng](){ return gsl_ran_flat(rng,1.,2.); }),
      //Locus 2: positions = beta(1,10) from 2 to 3
      std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,&pop.gametes[2],0.005,rng,[&rng](){ return 2. + gsl_ran_beta(rng,1.,10.); }),
      //Locus 2: positions = uniform [3,4)
      std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,&pop.gametes[3],0.005,rng,[&rng](){ return gsl_ran_flat(rng,3.,4.); })
      };

  //Equal mutation and rec. rates per locus
  std::vector<double> mu(4,0.005),rbw(3,0.005);
  BOOST_REQUIRE_EQUAL(mu.size(),4);
  for( ; generation < 10 ; ++generation )
    {
      double wbar = KTfwd::sample_diploid( rng,
					   &pop.gametes,
					   &pop.diploids,
					   &pop.mutations,
					   1000,
					   &mu[0],
					   mutmodels,
					   recmodels,
					   &rbw[0],
					   [](gsl_rng * __r, const double __d){ return gsl_ran_binomial(__r,__d,1); },
					   std::bind(KTfwd::insert_at_end<poptype::mutation_t,poptype::mlist_t>,std::placeholders::_1,std::placeholders::_2),
					   std::bind(KTfwd::insert_at_end<poptype::gamete_t,poptype::glist_t>,std::placeholders::_1,std::placeholders::_2),
					   std::bind(no_selection_multi(),std::placeholders::_1),
					   std::bind(KTfwd::mutation_remover(),std::placeholders::_1,0,2000),
					   0.);
      for( unsigned i = 0 ; i < 4 ; ++i )
	{
	  assert( check_sum(pop.gametes[i],2000) );
	}
      KTfwd::remove_fixed_lost(&pop.mutations,&pop.fixations,&pop.fixation_times,&pop.mut_lookup,generation,2000);
    }
  KTfwd::serialize s;
  s(pop,mwriter());
  poptype pop2(0,0);
  KTfwd::deserialize()(pop2,s,mreader());

  //Compare the mutations
  for( auto m1 = pop.mutations.begin(),m2 = pop2.mutations.begin() ; m1 != pop.mutations.end() ; ++m1,++m2 )
    {
      BOOST_CHECK_EQUAL( m1->pos, m2->pos );
      BOOST_CHECK_EQUAL( m1->n, m2->n );
    }

  //Compare the gametes
  for( auto gloc1 = pop.gametes.begin(), gloc2 = pop2.gametes.begin() ; gloc1 != pop.gametes.end() ; ++gloc1,++gloc2 )
    {
      for( auto g1 = gloc1->begin(),g2 = gloc2->begin() ; g1 != gloc1->end() ; ++g1,++g2 )
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
    }

  //Compare the diploids
  for( auto d1 = pop.diploids.begin(), d2 = pop2.diploids.begin() ; d1 != pop.diploids.end() ; ++d1,++d2 )
    {
      auto gloc1 = pop.gametes.begin();
      auto gloc2 = pop2.gametes.begin();
      //Iterate over loci w/in diploid
      for (auto l1 = d1->begin(), l2 = d2->begin() ; l1 != d1->end() ; ++l1,++l2)
	{
	  for( auto m1 = l1->first->mutations.begin(),m2=l2->first->mutations.begin() ; m1 != l1->first->mutations.end() ; ++m1,++m2 )
	    {
	      BOOST_CHECK( m1 != m2 );
	      BOOST_CHECK( std::distance(pop.mutations.begin(),*m1) == std::distance(pop2.mutations.begin(),*m2) );
	      BOOST_CHECK_EQUAL( (*m1)->pos, (*m2)->pos );
	      BOOST_CHECK_EQUAL( (*m1)->n, (*m2)->n ); 
	    }
	  for( auto m1 = l1->second->mutations.begin(),m2=l2->second->mutations.begin() ; m1 != l1->second->mutations.end() ; ++m1,++m2 )
	    {
	      BOOST_CHECK( m1 != m2 );
	      BOOST_CHECK( std::distance(pop.mutations.begin(),*m1) == std::distance(pop2.mutations.begin(),*m2) );
	      BOOST_CHECK_EQUAL( (*m1)->pos, (*m2)->pos );
	      BOOST_CHECK_EQUAL( (*m1)->n, (*m2)->n ); 
	    }
	}
    }
}

BOOST_AUTO_TEST_CASE( metapop_sugar_copy_construct )
{
  using poptype = KTfwd::multiloc_serialized<KTfwd::popgenmut,mwriter,mreader>;
  KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);

  poptype pop(1000,4);
  BOOST_REQUIRE_EQUAL( pop.gametes[0].size(),1 );
  BOOST_REQUIRE_EQUAL( pop.gametes[1].size(),1 );
  BOOST_REQUIRE_EQUAL( pop.gametes[2].size(),1 );
  BOOST_REQUIRE_EQUAL( pop.gametes[3].size(),1 );
  BOOST_REQUIRE_EQUAL( pop.diploids.size(), 1000 );
  for( unsigned i = 0 ; i < pop.diploids.size() ; ++i )
    {
      BOOST_REQUIRE_EQUAL( pop.diploids[i].size(), 4 );
      auto gptr = pop.gametes.begin();
      for( unsigned j = 0 ; j < 4 ; ++j, ++gptr )
	{
	  BOOST_REQUIRE( pop.diploids[i][j].first == gptr->begin() );
	  BOOST_REQUIRE( pop.diploids[i][j].second == gptr->begin() );
	}
    }
  BOOST_REQUIRE( pop.mutations.empty() == true );
  unsigned generation = 0;
  
  //Mutation model for 4 loci
  std::vector< std::function<KTfwd::popgenmut(typename poptype::mlist_t *)> > mutmodels {
    //Locus 0: positions Uniform [0,1)
    std::bind(KTfwd::infsites(),rng,std::placeholders::_1,&pop.mut_lookup,&generation,
	      0.005,0.,[](gsl_rng * r){return gsl_rng_uniform(r);},[](gsl_rng * r){return 0.;},[](gsl_rng * r){return 0.;}) ,
      //Locus 1: positions Uniform [1,2)
      std::bind(KTfwd::infsites(),rng,std::placeholders::_1,&pop.mut_lookup,&generation,
      		0.005,0.,[](gsl_rng * r){return gsl_ran_flat(r,1.,2.);},[](gsl_rng * r){return 0.;},[](gsl_rng * r){return 0.;}),
      //Locus 2: positions Uniform [2,3)
      std::bind(KTfwd::infsites(),rng,std::placeholders::_1,&pop.mut_lookup,&generation,
      		0.005,0.,[](gsl_rng * r){return gsl_ran_flat(r,2.,3.);},[](gsl_rng * r){return 0.;},[](gsl_rng * r){return 0.;}),
      //Locus 3: positions Uniform [3,4)
      std::bind(KTfwd::infsites(),rng,std::placeholders::_1,&pop.mut_lookup,&generation,
      		0.005,0.,[](gsl_rng * r){return gsl_ran_flat(r,3.,4.);},[](gsl_rng * r){return 0.;},[](gsl_rng * r){return 0.;})
  };
  
  //Within-locus recombination models for 4 loci
  std::vector< std::function<unsigned(typename poptype::glist_t::iterator &,
				      typename poptype::glist_t::iterator &)> > recmodels {
    //Locus 0: positions = uniform [0,1)
    std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,&pop.gametes[0],0.005,rng,[&rng](){ return gsl_rng_uniform(rng); }),
      //Locus 1: positions = uniform [1,2)
      std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,&pop.gametes[1],0.005,rng,[&rng](){ return gsl_ran_flat(rng,1.,2.); }),
      //Locus 2: positions = beta(1,10) from 2 to 3
      std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,&pop.gametes[2],0.005,rng,[&rng](){ return 2. + gsl_ran_beta(rng,1.,10.); }),
      //Locus 2: positions = uniform [3,4)
      std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,&pop.gametes[3],0.005,rng,[&rng](){ return gsl_ran_flat(rng,3.,4.); })
      };

  //Equal mutation and rec. rates per locus
  std::vector<double> mu(4,0.005),rbw(3,0.005);
  BOOST_REQUIRE_EQUAL(mu.size(),4);
  for( ; generation < 10 ; ++generation )
    {
      double wbar = KTfwd::sample_diploid( rng,
					   &pop.gametes,
					   &pop.diploids,
					   &pop.mutations,
					   1000,
					   &mu[0],
					   mutmodels,
					   recmodels,
					   &rbw[0],
					   [](gsl_rng * __r, const double __d){ return gsl_ran_binomial(__r,__d,1); },
					   std::bind(KTfwd::insert_at_end<poptype::mutation_t,poptype::mlist_t>,std::placeholders::_1,std::placeholders::_2),
					   std::bind(KTfwd::insert_at_end<poptype::gamete_t,poptype::glist_t>,std::placeholders::_1,std::placeholders::_2),
					   std::bind(no_selection_multi(),std::placeholders::_1),
					   std::bind(KTfwd::mutation_remover(),std::placeholders::_1,0,2000),
					   0.);
      for( unsigned i = 0 ; i < 4 ; ++i )
	{
	  assert( check_sum(pop.gametes[i],2000) );
	}
      KTfwd::remove_fixed_lost(&pop.mutations,&pop.fixations,&pop.fixation_times,&pop.mut_lookup,generation,2000);
    }

  poptype pop2(pop);

  //Compare the mutations
  for( auto m1 = pop.mutations.begin(),m2 = pop2.mutations.begin() ; m1 != pop.mutations.end() ; ++m1,++m2 )
    {
      BOOST_CHECK_EQUAL( m1->pos, m2->pos );
      BOOST_CHECK_EQUAL( m1->n, m2->n );
    }

  //Compare the gametes
  for( auto gloc1 = pop.gametes.begin(), gloc2 = pop2.gametes.begin() ; gloc1 != pop.gametes.end() ; ++gloc1,++gloc2 )
    {
      for( auto g1 = gloc1->begin(),g2 = gloc2->begin() ; g1 != gloc1->end() ; ++g1,++g2 )
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
    }

  //Compare the diploids
  for( auto d1 = pop.diploids.begin(), d2 = pop2.diploids.begin() ; d1 != pop.diploids.end() ; ++d1,++d2 )
    {
      auto gloc1 = pop.gametes.begin();
      auto gloc2 = pop2.gametes.begin();
      //Iterate over loci w/in diploid
      for (auto l1 = d1->begin(), l2 = d2->begin() ; l1 != d1->end() ; ++l1,++l2)
	{
	  for( auto m1 = l1->first->mutations.begin(),m2=l2->first->mutations.begin() ; m1 != l1->first->mutations.end() ; ++m1,++m2 )
	    {
	      BOOST_CHECK( m1 != m2 );
	      BOOST_CHECK( std::distance(pop.mutations.begin(),*m1) == std::distance(pop2.mutations.begin(),*m2) );
	      BOOST_CHECK_EQUAL( (*m1)->pos, (*m2)->pos );
	      BOOST_CHECK_EQUAL( (*m1)->n, (*m2)->n ); 
	    }
	  for( auto m1 = l1->second->mutations.begin(),m2=l2->second->mutations.begin() ; m1 != l1->second->mutations.end() ; ++m1,++m2 )
	    {
	      BOOST_CHECK( m1 != m2 );
	      BOOST_CHECK( std::distance(pop.mutations.begin(),*m1) == std::distance(pop2.mutations.begin(),*m2) );
	      BOOST_CHECK_EQUAL( (*m1)->pos, (*m2)->pos );
	      BOOST_CHECK_EQUAL( (*m1)->n, (*m2)->n ); 
	    }
	}
    }
}

BOOST_AUTO_TEST_CASE( metapop_sugar_assigment_operator )
{
  using poptype = KTfwd::multiloc_serialized<KTfwd::popgenmut,mwriter,mreader>;
  KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);

  poptype pop(1000,4);
  BOOST_REQUIRE_EQUAL( pop.gametes[0].size(),1 );
  BOOST_REQUIRE_EQUAL( pop.gametes[1].size(),1 );
  BOOST_REQUIRE_EQUAL( pop.gametes[2].size(),1 );
  BOOST_REQUIRE_EQUAL( pop.gametes[3].size(),1 );
  BOOST_REQUIRE_EQUAL( pop.diploids.size(), 1000 );
  for( unsigned i = 0 ; i < pop.diploids.size() ; ++i )
    {
      BOOST_REQUIRE_EQUAL( pop.diploids[i].size(), 4 );
      auto gptr = pop.gametes.begin();
      for( unsigned j = 0 ; j < 4 ; ++j, ++gptr )
	{
	  BOOST_REQUIRE( pop.diploids[i][j].first == gptr->begin() );
	  BOOST_REQUIRE( pop.diploids[i][j].second == gptr->begin() );
	}
    }
  BOOST_REQUIRE( pop.mutations.empty() == true );
  unsigned generation = 0;
  
  //Mutation model for 4 loci
  std::vector< std::function<KTfwd::popgenmut(typename poptype::mlist_t *)> > mutmodels {
    //Locus 0: positions Uniform [0,1)
    std::bind(KTfwd::infsites(),rng,std::placeholders::_1,&pop.mut_lookup,&generation,
	      0.005,0.,[](gsl_rng * r){return gsl_rng_uniform(r);},[](gsl_rng * r){return 0.;},[](gsl_rng * r){return 0.;}) ,
      //Locus 1: positions Uniform [1,2)
      std::bind(KTfwd::infsites(),rng,std::placeholders::_1,&pop.mut_lookup,&generation,
      		0.005,0.,[](gsl_rng * r){return gsl_ran_flat(r,1.,2.);},[](gsl_rng * r){return 0.;},[](gsl_rng * r){return 0.;}),
      //Locus 2: positions Uniform [2,3)
      std::bind(KTfwd::infsites(),rng,std::placeholders::_1,&pop.mut_lookup,&generation,
      		0.005,0.,[](gsl_rng * r){return gsl_ran_flat(r,2.,3.);},[](gsl_rng * r){return 0.;},[](gsl_rng * r){return 0.;}),
      //Locus 3: positions Uniform [3,4)
      std::bind(KTfwd::infsites(),rng,std::placeholders::_1,&pop.mut_lookup,&generation,
      		0.005,0.,[](gsl_rng * r){return gsl_ran_flat(r,3.,4.);},[](gsl_rng * r){return 0.;},[](gsl_rng * r){return 0.;})
  };
  
  //Within-locus recombination models for 4 loci
  std::vector< std::function<unsigned(typename poptype::glist_t::iterator &,
				      typename poptype::glist_t::iterator &)> > recmodels {
    //Locus 0: positions = uniform [0,1)
    std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,&pop.gametes[0],0.005,rng,[&rng](){ return gsl_rng_uniform(rng); }),
      //Locus 1: positions = uniform [1,2)
      std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,&pop.gametes[1],0.005,rng,[&rng](){ return gsl_ran_flat(rng,1.,2.); }),
      //Locus 2: positions = beta(1,10) from 2 to 3
      std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,&pop.gametes[2],0.005,rng,[&rng](){ return 2. + gsl_ran_beta(rng,1.,10.); }),
      //Locus 2: positions = uniform [3,4)
      std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,&pop.gametes[3],0.005,rng,[&rng](){ return gsl_ran_flat(rng,3.,4.); })
      };

  //Equal mutation within and rec. rates b/w loci
  std::vector<double> mu(4,0.005),rbw(3,0.005);
  BOOST_REQUIRE_EQUAL(mu.size(),4);
  for( ; generation < 10 ; ++generation )
    {
      double wbar = KTfwd::sample_diploid( rng,
					   &pop.gametes,
					   &pop.diploids,
					   &pop.mutations,
					   1000,
					   &mu[0],
					   mutmodels,
					   recmodels,
					   &rbw[0],
					   [](gsl_rng * __r, const double __d){ return gsl_ran_binomial(__r,__d,1); },
					   std::bind(KTfwd::insert_at_end<poptype::mutation_t,poptype::mlist_t>,std::placeholders::_1,std::placeholders::_2),
					   std::bind(KTfwd::insert_at_end<poptype::gamete_t,poptype::glist_t>,std::placeholders::_1,std::placeholders::_2),
					   std::bind(no_selection_multi(),std::placeholders::_1),
					   std::bind(KTfwd::mutation_remover(),std::placeholders::_1,0,2000),
					   0.);
      for( unsigned i = 0 ; i < 4 ; ++i )
	{
	  assert( check_sum(pop.gametes[i],2000) );
	}
      KTfwd::remove_fixed_lost(&pop.mutations,&pop.fixations,&pop.fixation_times,&pop.mut_lookup,generation,2000);
    }

  poptype pop2 = pop;

  //Compare the mutations
  for( auto m1 = pop.mutations.begin(),m2 = pop2.mutations.begin() ; m1 != pop.mutations.end() ; ++m1,++m2 )
    {
      BOOST_CHECK_EQUAL( m1->pos, m2->pos );
      BOOST_CHECK_EQUAL( m1->n, m2->n );
    }

  //Compare the gametes
  for( auto gloc1 = pop.gametes.begin(), gloc2 = pop2.gametes.begin() ; gloc1 != pop.gametes.end() ; ++gloc1,++gloc2 )
    {
      for( auto g1 = gloc1->begin(),g2 = gloc2->begin() ; g1 != gloc1->end() ; ++g1,++g2 )
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
    }

  //Compare the diploids
  for( auto d1 = pop.diploids.begin(), d2 = pop2.diploids.begin() ; d1 != pop.diploids.end() ; ++d1,++d2 )
    {
      auto gloc1 = pop.gametes.begin();
      auto gloc2 = pop2.gametes.begin();
      //Iterate over loci w/in diploid
      for (auto l1 = d1->begin(), l2 = d2->begin() ; l1 != d1->end() ; ++l1,++l2)
	{
	  for( auto m1 = l1->first->mutations.begin(),m2=l2->first->mutations.begin() ; m1 != l1->first->mutations.end() ; ++m1,++m2 )
	    {
	      BOOST_CHECK( m1 != m2 );
	      BOOST_CHECK( std::distance(pop.mutations.begin(),*m1) == std::distance(pop2.mutations.begin(),*m2) );
	      BOOST_CHECK_EQUAL( (*m1)->pos, (*m2)->pos );
	      BOOST_CHECK_EQUAL( (*m1)->n, (*m2)->n ); 
	    }
	  for( auto m1 = l1->second->mutations.begin(),m2=l2->second->mutations.begin() ; m1 != l1->second->mutations.end() ; ++m1,++m2 )
	    {
	      BOOST_CHECK( m1 != m2 );
	      BOOST_CHECK( std::distance(pop.mutations.begin(),*m1) == std::distance(pop2.mutations.begin(),*m2) );
	      BOOST_CHECK_EQUAL( (*m1)->pos, (*m2)->pos );
	      BOOST_CHECK_EQUAL( (*m1)->n, (*m2)->n ); 
	    }
	}
    }
}
