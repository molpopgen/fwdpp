/*!
  \file sugar3.cc
  \ingroup unit
  \brief Testing KTfwd::multiloc and KTfwd::multiloc_serialized
*/

#define BOOST_TEST_MODULE sugarTest3
#define BOOST_TEST_DYN_LINK 

#include <iostream>
#include <unistd.h>
#include <config.h>
#include <functional>
#include <algorithm>
#include <boost/test/unit_test.hpp>
#include <fwdpp/diploid.hh>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/sugar/multiloc.hpp>
#include <fwdpp/sugar/infsites.hpp>

using mwriter = KTfwd::mutation_writer;
using mreader = KTfwd::mutation_reader<KTfwd::popgenmut>;

//Fitness function
struct no_selection_multi
{
  using poptype = KTfwd::multiloc<KTfwd::popgenmut>;
  using result_type = double;
  inline double operator()(const poptype::dipvector_t::const_iterator ) const
  {
    return 1.;
  }
};

BOOST_AUTO_TEST_CASE( multiloc_sugar_test1 )
{
  using poptype = KTfwd::multiloc<KTfwd::popgenmut>;
  KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);

  poptype pop(1000,4);
  BOOST_REQUIRE_EQUAL( pop.gametes.size(),1 );
  BOOST_REQUIRE_EQUAL( pop.diploids.size(), 1000 );
  for( unsigned i = 0 ; i < pop.diploids.size() ; ++i )
    {
      BOOST_REQUIRE_EQUAL( pop.diploids[i].size(), 4 );
      auto gptr = pop.gametes.begin();
      for( unsigned j = 0 ; j < 4 ; ++j )
	{
	  BOOST_REQUIRE( pop.diploids[i][j].first == gptr );
	  BOOST_REQUIRE( pop.diploids[i][j].second == gptr );
	}
    }
  BOOST_REQUIRE( pop.mutations.empty() == true );
  unsigned generation = 0;
  
  //Mutation model for 4 loci
  //auto x = std::result_of<KTfwd::fwdpp_internal::make_mut_queue(poptype::mlist_t *)>::type();
  
  std::vector< std::function<poptype::mlist_t::iterator(decltype(KTfwd::fwdpp_internal::make_mut_queue(&pop.mutations)) &,
							typename poptype::mlist_t *)> > mutmodels {
    //Locus 0: positions Uniform [0,1)
    std::bind(KTfwd::infsites(),std::placeholders::_1,std::placeholders::_2,rng.get(),&pop.mut_lookup,&generation,
	      0.005,0.,[&rng](){return gsl_rng_uniform(rng.get());},[](){return 0.;},[](){return 0.;}) ,
      //Locus 1: positions Uniform [1,2)
      std::bind(KTfwd::infsites(),std::placeholders::_1,std::placeholders::_2,rng.get(),&pop.mut_lookup,&generation,
      		0.005,0.,[&rng](){return gsl_ran_flat(rng.get(),1.,2.);},[](){return 0.;},[](){return 0.;}),
      //Locus 2: positions Uniform [2,3)
      std::bind(KTfwd::infsites(),std::placeholders::_1,std::placeholders::_2,rng.get(),&pop.mut_lookup,&generation,
      		0.005,0.,[&rng](){return gsl_ran_flat(rng.get(),2.,3.);},[](){return 0.;},[](){return 0.;}),
      //Locus 3: positions Uniform [3,4)
      std::bind(KTfwd::infsites(),std::placeholders::_1,std::placeholders::_2,rng.get(),&pop.mut_lookup,&generation,
      		0.005,0.,[&rng](){return gsl_ran_flat(rng.get(),3.,4.);},[](){return 0.;},[](){return 0.;})
      };
  
  //Within-locus recombination models for 4 loci
  std::vector< std::function<unsigned(typename poptype::glist_t::iterator &,
				      typename poptype::glist_t::iterator &,
				      decltype(KTfwd::fwdpp_internal::gamete_lookup_table(&pop.gametes)) &,
				      decltype(KTfwd::fwdpp_internal::make_gamete_queue(&pop.gametes)) &)
			     > > recmodels {
    //Locus 0: positions = uniform [0,1)
    std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,std::placeholders::_4,std::ref(pop.neutral),std::ref(pop.selected),&pop.gametes,0.005,rng.get(),[&rng](){ return gsl_rng_uniform(rng.get()); }),
      //Locus 1: positions = uniform [1,2)
      std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,std::placeholders::_4,std::ref(pop.neutral),std::ref(pop.selected),&pop.gametes,0.005,rng.get(),[&rng](){ return gsl_ran_flat(rng.get(),1.,2.); }),
      //Locus 2: positions = beta(1,10) from 2 to 3
      std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,std::placeholders::_4,std::ref(pop.neutral),std::ref(pop.selected),&pop.gametes,0.005,rng.get(),[&rng](){ return 2. + gsl_ran_beta(rng.get(),1.,10.); }),
      //Locus 2: positions = uniform [3,4)
      std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,std::placeholders::_4,std::ref(pop.neutral),std::ref(pop.selected),&pop.gametes,0.005,rng.get(),[&rng](){ return gsl_ran_flat(rng.get(),3.,4.); })
      };

  //Equal mutation and rec. rates per locus
  std::vector<double> mu(4,0.005),rbw(3,0.005);
  BOOST_REQUIRE_EQUAL(mu.size(),4);
  for( ; generation < 10 ; ++generation )
    {
      double wbar = KTfwd::sample_diploid( rng.get(),
					   &pop.gametes,
					   &pop.diploids,
					   &pop.mutations,
					   1000,
					   &mu[0],
					   mutmodels,
					   recmodels,
					   &rbw[0],
					   [](gsl_rng * __r, const double __d){ return gsl_ran_binomial(__r,__d,1); },
					   std::bind(KTfwd::insert_at_end<poptype::gamete_t,poptype::glist_t>,std::placeholders::_1,std::placeholders::_2),
					   std::bind(no_selection_multi(),std::placeholders::_1),
					   std::bind(KTfwd::mutation_remover(),std::placeholders::_1,2000),
					   0.);
      assert( check_sum(pop.gametes,8000) );
      KTfwd::update_mutations(&pop.mutations,&pop.fixations,&pop.fixation_times,&pop.mut_lookup,generation,2000);
    }
  KTfwd::serialize s;
  s(pop,mwriter());
  poptype pop2(0,0);
  KTfwd::deserialize()(pop2,s,mreader());

  BOOST_REQUIRE_EQUAL( pop.mutations.size(),pop2.mutations.size() );
  BOOST_REQUIRE_EQUAL( pop.gametes.size(),pop2.gametes.size() );
  BOOST_REQUIRE_EQUAL( pop.diploids.size(),pop2.diploids.size() );
  for( unsigned i=0;i<pop.diploids.size();++i) BOOST_REQUIRE_EQUAL( pop.diploids[i].size(),pop2.diploids[i].size() );
  //Compare the mutations
  for( auto m1 = pop.mutations.begin(),m2 = pop2.mutations.begin() ; m1 != pop.mutations.end() ; ++m1,++m2 )
    {
      BOOST_CHECK_EQUAL( m1->pos, m2->pos );
      BOOST_CHECK_EQUAL( m1->n, m2->n );
    }

  //Compare the gametes
  for( auto gloc1 = pop.gametes.begin(), gloc2 = pop2.gametes.begin() ; gloc1 != pop.gametes.end() ; ++gloc1,++gloc2 )
    {
      for( auto g1 = pop.gametes.begin(),g2 = pop2.gametes.begin() ; g1 != pop.gametes.end() ; ++g1,++g2 )
  	{
  	  BOOST_CHECK( g1 != g2 );
	  if(g1->n) {
	    for( auto m1 = g1->mutations.begin(),m2=g2->mutations.begin() ; m1 != g1->mutations.end() ; ++m1,++m2 )
	      {
		BOOST_CHECK( m1 != m2 );
		BOOST_CHECK( std::distance(pop.mutations.begin(),*m1) == std::distance(pop2.mutations.begin(),*m2) );
		BOOST_CHECK_EQUAL( (*m1)->pos, (*m2)->pos );
		BOOST_CHECK_EQUAL( (*m1)->n, (*m2)->n ); 
	      }
	  }
  	}
    }

  //Compare the diploids
  for( auto d1 = pop.diploids.begin(), d2 = pop2.diploids.begin() ; d1 != pop.diploids.end() ; ++d1,++d2 )
    {
      //Iterate over loci w/in diploid
      for (auto l1 = d1->begin(), l2 = d2->begin() ; l1 != d1->end() ; ++l1,++l2)
	{
	  BOOST_REQUIRE_EQUAL(l1->first->mutations.size(),l2->first->mutations.size());
	  BOOST_REQUIRE_EQUAL(l1->second->mutations.size(),l2->second->mutations.size());
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

BOOST_AUTO_TEST_CASE( multiloc_sugar_gzserialize_test )
{
  using poptype = KTfwd::multiloc<KTfwd::popgenmut>;
  KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);

  poptype pop(1000,4);
  BOOST_REQUIRE_EQUAL( pop.gametes.size(),1 );
  BOOST_REQUIRE_EQUAL( pop.diploids.size(), 1000 );
  for( unsigned i = 0 ; i < pop.diploids.size() ; ++i )
    {
      BOOST_REQUIRE_EQUAL( pop.diploids[i].size(), 4 );
      auto gptr = pop.gametes.begin();
      for( unsigned j = 0 ; j < 4 ; ++j )
	{
	  BOOST_REQUIRE( pop.diploids[i][j].first == gptr );
	  BOOST_REQUIRE( pop.diploids[i][j].second == gptr );
	}
    }
  BOOST_REQUIRE( pop.mutations.empty() == true );
  unsigned generation = 0;
  
  //Mutation model for 4 loci
  std::vector< std::function<poptype::mlist_t::iterator(decltype(KTfwd::fwdpp_internal::make_mut_queue(&pop.mutations)) &,
							typename poptype::mlist_t *)> > mutmodels {
    //Locus 0: positions Uniform [0,1)
    std::bind(KTfwd::infsites(),std::placeholders::_1,std::placeholders::_2,rng.get(),&pop.mut_lookup,&generation,
	      0.005,0.,[&rng](){return gsl_rng_uniform(rng.get());},[](){return 0.;},[](){return 0.;}) ,
      //Locus 1: positions Uniform [1,2)
      std::bind(KTfwd::infsites(),std::placeholders::_1,std::placeholders::_2,rng.get(),&pop.mut_lookup,&generation,
      		0.005,0.,[&rng](){return gsl_ran_flat(rng.get(),1.,2.);},[](){return 0.;},[](){return 0.;}),
      //Locus 2: positions Uniform [2,3)
      std::bind(KTfwd::infsites(),std::placeholders::_1,std::placeholders::_2,rng.get(),&pop.mut_lookup,&generation,
      		0.005,0.,[&rng](){return gsl_ran_flat(rng.get(),2.,3.);},[](){return 0.;},[](){return 0.;}),
      //Locus 3: positions Uniform [3,4)
      std::bind(KTfwd::infsites(),std::placeholders::_1,std::placeholders::_2,rng.get(),&pop.mut_lookup,&generation,
      		0.005,0.,[&rng](){return gsl_ran_flat(rng.get(),3.,4.);},[](){return 0.;},[](){return 0.;})
      };
  
  //Within-locus recombination models for 4 loci
  std::vector< std::function<unsigned(typename poptype::glist_t::iterator &,
				      typename poptype::glist_t::iterator &,
				      decltype(KTfwd::fwdpp_internal::gamete_lookup_table(&pop.gametes)) &,
				      decltype(KTfwd::fwdpp_internal::make_gamete_queue(&pop.gametes)) &)
			     > > recmodels {
    //Locus 0: positions = uniform [0,1)
    std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,std::placeholders::_4,std::ref(pop.neutral),std::ref(pop.selected),&pop.gametes,0.005,rng.get(),[&rng](){ return gsl_rng_uniform(rng.get()); }),
      //Locus 1: positions = uniform [1,2)
      std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,std::placeholders::_4,std::ref(pop.neutral),std::ref(pop.selected),&pop.gametes,0.005,rng.get(),[&rng](){ return gsl_ran_flat(rng.get(),1.,2.); }),
      //Locus 2: positions = beta(1,10) from 2 to 3
      std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,std::placeholders::_4,std::ref(pop.neutral),std::ref(pop.selected),&pop.gametes,0.005,rng.get(),[&rng](){ return 2. + gsl_ran_beta(rng.get(),1.,10.); }),
      //Locus 2: positions = uniform [3,4)
      std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,std::placeholders::_4,std::ref(pop.neutral),std::ref(pop.selected),&pop.gametes,0.005,rng.get(),[&rng](){ return gsl_ran_flat(rng.get(),3.,4.); })
      };

  //Equal mutation and rec. rates per locus
  std::vector<double> mu(4,0.005),rbw(3,0.005);
  BOOST_REQUIRE_EQUAL(mu.size(),4);
  for( ; generation < 10 ; ++generation )
    {
      double wbar = KTfwd::sample_diploid( rng.get(),
					   &pop.gametes,
					   &pop.diploids,
					   &pop.mutations,
					   1000,
					   &mu[0],
					   mutmodels,
					   recmodels,
					   &rbw[0],
					   [](gsl_rng * __r, const double __d){ return gsl_ran_binomial(__r,__d,1); },
					   std::bind(KTfwd::insert_at_end<poptype::gamete_t,poptype::glist_t>,std::placeholders::_1,std::placeholders::_2),
					   std::bind(no_selection_multi(),std::placeholders::_1),
					   std::bind(KTfwd::mutation_remover(),std::placeholders::_1,2000),
					   0.);
      assert( check_sum(pop.gametes,8000) );
      KTfwd::update_mutations(&pop.mutations,&pop.fixations,&pop.fixation_times,&pop.mut_lookup,generation,2000);
    }
  KTfwd::gzserialize s;
  gzFile gzf = gzopen("sugar3_test.gz","wb");
  s(gzf,pop,mwriter());
  gzclose(gzf);
  gzf = gzopen("sugar3_test.gz","rb");
  poptype pop2(0,0);
  KTfwd::gzdeserialize()(gzf,pop2,mreader());

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
      if(g1->n) {
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
  unlink("sugar3_test.gz");
}
BOOST_AUTO_TEST_CASE( multiloc_sugar_copy_construct )
{
  using poptype = KTfwd::multiloc_serialized<KTfwd::popgenmut,mwriter,mreader>;
  KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);

  poptype pop(1000,4);
  BOOST_REQUIRE_EQUAL( pop.gametes.size(),1 );
  BOOST_REQUIRE_EQUAL( pop.diploids.size(), 1000 );
  for( unsigned i = 0 ; i < pop.diploids.size() ; ++i )
    {
      BOOST_REQUIRE_EQUAL( pop.diploids[i].size(), 4 );
      auto gptr = pop.gametes.begin();
      for( unsigned j = 0 ; j < 4 ; ++j )
	{
	  BOOST_REQUIRE( pop.diploids[i][j].first == gptr );
	  BOOST_REQUIRE( pop.diploids[i][j].second == gptr );
	}
    }
  BOOST_REQUIRE( pop.mutations.empty() == true );
  unsigned generation = 0;
  
  //Mutation model for 4 loci
  std::vector< std::function<poptype::mlist_t::iterator(decltype(KTfwd::fwdpp_internal::make_mut_queue(&pop.mutations)) &,
							typename poptype::mlist_t *)> > mutmodels {
    //Locus 0: positions Uniform [0,1)
    std::bind(KTfwd::infsites(),std::placeholders::_1,std::placeholders::_2,rng.get(),&pop.mut_lookup,&generation,
	      0.005,0.,[&rng](){return gsl_rng_uniform(rng.get());},[](){return 0.;},[](){return 0.;}) ,
      //Locus 1: positions Uniform [1,2)
      std::bind(KTfwd::infsites(),std::placeholders::_1,std::placeholders::_2,rng.get(),&pop.mut_lookup,&generation,
      		0.005,0.,[&rng](){return gsl_ran_flat(rng.get(),1.,2.);},[](){return 0.;},[](){return 0.;}),
      //Locus 2: positions Uniform [2,3)
      std::bind(KTfwd::infsites(),std::placeholders::_1,std::placeholders::_2,rng.get(),&pop.mut_lookup,&generation,
      		0.005,0.,[&rng](){return gsl_ran_flat(rng.get(),2.,3.);},[](){return 0.;},[](){return 0.;}),
      //Locus 3: positions Uniform [3,4)
      std::bind(KTfwd::infsites(),std::placeholders::_1,std::placeholders::_2,rng.get(),&pop.mut_lookup,&generation,
      		0.005,0.,[&rng](){return gsl_ran_flat(rng.get(),3.,4.);},[](){return 0.;},[](){return 0.;})
      };
  
  //Within-locus recombination models for 4 loci
  std::vector< std::function<unsigned(typename poptype::glist_t::iterator &,
				      typename poptype::glist_t::iterator &,
				      decltype(KTfwd::fwdpp_internal::gamete_lookup_table(&pop.gametes)) &,
				      decltype(KTfwd::fwdpp_internal::make_gamete_queue(&pop.gametes)) &)
			     > > recmodels {
    //Locus 0: positions = uniform [0,1)
    std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,std::placeholders::_4,std::ref(pop.neutral),std::ref(pop.selected),&pop.gametes,0.005,rng.get(),[&rng](){ return gsl_rng_uniform(rng.get()); }),
      //Locus 1: positions = uniform [1,2)
      std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,std::placeholders::_4,std::ref(pop.neutral),std::ref(pop.selected),&pop.gametes,0.005,rng.get(),[&rng](){ return gsl_ran_flat(rng.get(),1.,2.); }),
      //Locus 2: positions = beta(1,10) from 2 to 3
      std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,std::placeholders::_4,std::ref(pop.neutral),std::ref(pop.selected),&pop.gametes,0.005,rng.get(),[&rng](){ return 2. + gsl_ran_beta(rng.get(),1.,10.); }),
      //Locus 2: positions = uniform [3,4)
      std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,std::placeholders::_4,std::ref(pop.neutral),std::ref(pop.selected),&pop.gametes,0.005,rng.get(),[&rng](){ return gsl_ran_flat(rng.get(),3.,4.); })
      };

  //Equal mutation and rec. rates per locus
  std::vector<double> mu(4,0.005),rbw(3,0.005);
  BOOST_REQUIRE_EQUAL(mu.size(),4);
  for( ; generation < 10 ; ++generation )
    {
      double wbar = KTfwd::sample_diploid( rng.get(),
					   &pop.gametes,
					   &pop.diploids,
					   &pop.mutations,
					   1000,
					   &mu[0],
					   mutmodels,
					   recmodels,
					   &rbw[0],
					   [](gsl_rng * __r, const double __d){ return gsl_ran_binomial(__r,__d,1); },
					   std::bind(KTfwd::insert_at_end<poptype::gamete_t,poptype::glist_t>,std::placeholders::_1,std::placeholders::_2),
					   std::bind(no_selection_multi(),std::placeholders::_1),
					   std::bind(KTfwd::mutation_remover(),std::placeholders::_1,2000),
					   0.);
      assert( check_sum(pop.gametes,8000) );
      KTfwd::update_mutations(&pop.mutations,&pop.fixations,&pop.fixation_times,&pop.mut_lookup,generation,2000);
    }

  poptype pop2(pop);

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
      if(g1->n){
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

BOOST_AUTO_TEST_CASE( multiloc_sugar_assigment_operator )
{
  using poptype = KTfwd::multiloc_serialized<KTfwd::popgenmut,mwriter,mreader>;
  KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);

  poptype pop(1000,4);
  BOOST_REQUIRE_EQUAL( pop.gametes.size(),1 );
  BOOST_REQUIRE_EQUAL( pop.diploids.size(), 1000 );
  for( unsigned i = 0 ; i < pop.diploids.size() ; ++i )
    {
      BOOST_REQUIRE_EQUAL( pop.diploids[i].size(), 4 );
      auto gptr = pop.gametes.begin();
      for( unsigned j = 0 ; j < 4 ; ++j )
	{
	  BOOST_REQUIRE( pop.diploids[i][j].first == gptr);
	  BOOST_REQUIRE( pop.diploids[i][j].second == gptr );
	}
    }
  BOOST_REQUIRE( pop.mutations.empty() == true );
  unsigned generation = 0;
  
  //Mutation model for 4 loci
  std::vector< std::function<poptype::mlist_t::iterator(decltype(KTfwd::fwdpp_internal::make_mut_queue(&pop.mutations)) &,
							typename poptype::mlist_t *)> > mutmodels {
    //Locus 0: positions Uniform [0,1)
    std::bind(KTfwd::infsites(),std::placeholders::_1,std::placeholders::_2,rng.get(),&pop.mut_lookup,&generation,
	      0.005,0.,[&rng](){return gsl_rng_uniform(rng.get());},[](){return 0.;},[](){return 0.;}) ,
      //Locus 1: positions Uniform [1,2)
      std::bind(KTfwd::infsites(),std::placeholders::_1,std::placeholders::_2,rng.get(),&pop.mut_lookup,&generation,
      		0.005,0.,[&rng](){return gsl_ran_flat(rng.get(),1.,2.);},[](){return 0.;},[](){return 0.;}),
      //Locus 2: positions Uniform [2,3)
      std::bind(KTfwd::infsites(),std::placeholders::_1,std::placeholders::_2,rng.get(),&pop.mut_lookup,&generation,
      		0.005,0.,[&rng](){return gsl_ran_flat(rng.get(),2.,3.);},[](){return 0.;},[](){return 0.;}),
      //Locus 3: positions Uniform [3,4)
      std::bind(KTfwd::infsites(),std::placeholders::_1,std::placeholders::_2,rng.get(),&pop.mut_lookup,&generation,
      		0.005,0.,[&rng](){return gsl_ran_flat(rng.get(),3.,4.);},[](){return 0.;},[](){return 0.;})
      };
  
  //Within-locus recombination models for 4 loci
  std::vector< std::function<unsigned(typename poptype::glist_t::iterator &,
				      typename poptype::glist_t::iterator &,
				      decltype(KTfwd::fwdpp_internal::gamete_lookup_table(&pop.gametes)) &,
				      decltype(KTfwd::fwdpp_internal::make_gamete_queue(&pop.gametes)) &)> > recmodels {
    //Locus 0: positions = uniform [0,1)
    std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,std::placeholders::_4,std::ref(pop.neutral),std::ref(pop.selected),&pop.gametes,0.005,rng.get(),[&rng](){ return gsl_rng_uniform(rng.get()); }),
      //Locus 1: positions = uniform [1,2)
      std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,std::placeholders::_4,std::ref(pop.neutral),std::ref(pop.selected),&pop.gametes,0.005,rng.get(),[&rng](){ return gsl_ran_flat(rng.get(),1.,2.); }),
      //Locus 2: positions = beta(1,10) from 2 to 3
      std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,std::placeholders::_4,std::ref(pop.neutral),std::ref(pop.selected),&pop.gametes,0.005,rng.get(),[&rng](){ return 2. + gsl_ran_beta(rng.get(),1.,10.); }),
      //Locus 2: positions = uniform [3,4)
      std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,std::placeholders::_4,std::ref(pop.neutral),std::ref(pop.selected),&pop.gametes,0.005,rng.get(),[&rng](){ return gsl_ran_flat(rng.get(),3.,4.); })
      };

  //Equal mutation within and rec. rates b/w loci
  std::vector<double> mu(4,0.005),rbw(3,0.005);
  BOOST_REQUIRE_EQUAL(mu.size(),4);
  for( ; generation < 10 ; ++generation )
    {
      double wbar = KTfwd::sample_diploid( rng.get(),
					   &pop.gametes,
					   &pop.diploids,
					   &pop.mutations,
					   1000,
					   &mu[0],
					   mutmodels,
					   recmodels,
					   &rbw[0],
					   [](gsl_rng * __r, const double __d){ return gsl_ran_binomial(__r,__d,1); },
					   std::bind(KTfwd::insert_at_end<poptype::gamete_t,poptype::glist_t>,std::placeholders::_1,std::placeholders::_2),
					   std::bind(no_selection_multi(),std::placeholders::_1),
					   std::bind(KTfwd::mutation_remover(),std::placeholders::_1,2000),
					   0.);
      assert( check_sum(pop.gametes,8000) );
      KTfwd::update_mutations(&pop.mutations,&pop.fixations,&pop.fixation_times,&pop.mut_lookup,generation,2000);
    }

  poptype pop2 = pop;

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
      if(g1->n)
	{
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
