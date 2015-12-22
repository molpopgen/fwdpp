/*! 
  \file sugar1_custom.cc 
  \ingroup unit 
  \brief Testing single-deme sugar functionality with custom diploids
*/
#define BOOST_TEST_MODULE sugarTest1_custom
#define BOOST_TEST_DYN_LINK 

#include <config.h>
#include <boost/test/unit_test.hpp>
#include <fwdpp/diploid.hh>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/sugar/singlepop.hpp>
#include <fwdpp/sugar/infsites.hpp>

using mutation_t = KTfwd::popgenmut;
using mwriter = KTfwd::mutation_writer;
using mreader = KTfwd::mutation_reader<mutation_t>;

//Custom diploid type.
struct diploid_t : public KTfwd::tags::custom_diploid_t
{
  using first_type = KTfwd::singlepop_glist_t<mutation_t>::iterator;
  using second_type = KTfwd::singlepop_glist_t<mutation_t>::iterator;
  first_type first;
  second_type second;
  unsigned i;
  diploid_t() : first(first_type()),second(second_type()),i(std::numeric_limits<unsigned>::max()) {}
  diploid_t(first_type g1, first_type g2) : first(g1),second(g2),i(std::numeric_limits<unsigned>::max()) {}
};

struct diploid_writer
{
  using result_type = void;
  template<typename itr,typename streamtype>
  inline result_type operator()( itr i, streamtype & o ) const
  {
    o.write( reinterpret_cast<const char *>(&i->i),sizeof(unsigned) );
  }
};

struct diploid_reader
{
  using result_type = void;
  template<typename itr,typename streamtype>
  inline result_type operator()( itr i, streamtype & in ) const
  {
    in.read( reinterpret_cast<char *>(&i->i),sizeof(unsigned) );
  }
};

BOOST_AUTO_TEST_CASE( singlepop_sugar_test1 )
{
  using poptype = KTfwd::singlepop<mutation_t,diploid_t>;
  poptype pop(1000);

  KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);

  //Evolve for 10 generations
  std::function<double(void)> recmap = std::bind(gsl_rng_uniform,rng.get());
  for( unsigned generation= 0 ; generation < 10 ; ++generation )
    {
      double wbar = KTfwd::sample_diploid(rng.get(),
					  &pop.gametes,
					  &pop.diploids,
					  &pop.mutations,
					  1000,
					  0.005,
					  std::bind(KTfwd::infsites(),std::placeholders::_1,std::placeholders::_2,rng.get(),&pop.mut_lookup,generation,
						    0.005,0.,[&rng](){return gsl_rng_uniform(rng.get());},[](){return 0.;},[](){return 0.;}),
					  std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,std::placeholders::_4,
						    std::ref(pop.neutral),std::ref(pop.selected),
						    &pop.gametes,
						    0., //no rec
						    rng.get(),
						    recmap),
					  std::bind(KTfwd::insert_at_end<poptype::gamete_t,poptype::glist_t>,std::placeholders::_1,std::placeholders::_2),
					  []( poptype::dipvector_t::const_iterator dip ) { return KTfwd::multiplicative_diploid()(dip->first,dip->second,2.); },
					  std::bind(KTfwd::mutation_remover(),std::placeholders::_1,2*pop.N));
      KTfwd::update_mutations(&pop.mutations,&pop.fixations,&pop.fixation_times,&pop.mut_lookup,generation,2*pop.N);
    }

  //This is the test--we will add "ID" numbers to each diploid at the end of the sim
  //This has zero relevance to any model, and is just a test.
  unsigned i=0;
  std::for_each(pop.diploids.begin(),pop.diploids.end(),
		[&i]( diploid_t & dip ) { dip.i = i++; } );

  //attempt to copy
  KTfwd::serialize s;
  s(pop,mwriter(),diploid_writer());
  poptype pop2(0);
  KTfwd::deserialize()(pop2,s,mreader(),diploid_reader());
  BOOST_REQUIRE(pop.mutations.size() == pop2.mutations.size());
  // unsigned NN=0;
  // for(const auto & g : pop.gametes)
  //   {
  //     if(g.n)++NN;
  //   }
  // BOOST_REQUIRE(NN == pop2.gametes.size());
  BOOST_REQUIRE(pop.gametes.size() == pop2.gametes.size());
  BOOST_REQUIRE(pop.diploids.size() == pop2.diploids.size());
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
  i=0;
  for( auto d1 = pop.diploids.begin(),d2 = pop2.diploids.begin() ; d1 != pop.diploids.end() ; ++d1,++d2,++i )
    {
      BOOST_CHECK(d1->i == d2->i);//This is the serialized "extra" data!!
      BOOST_CHECK(d2->i==i);
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

BOOST_AUTO_TEST_CASE( singlepop_sugar_serialize_in_memory )
{
  using poptype = KTfwd::singlepop<mutation_t,diploid_t>;
  poptype pop(1000);

  KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);

  //Evolve for 10 generations
  std::function<double(void)> recmap = std::bind(gsl_rng_uniform,rng.get());
  for( unsigned generation= 0 ; generation < 10 ; ++generation )
    {
      double wbar = KTfwd::sample_diploid(rng.get(),
					  &pop.gametes,
					  &pop.diploids,
					  &pop.mutations,
					  1000,
					  0.005,
					  std::bind(KTfwd::infsites(),std::placeholders::_1,std::placeholders::_2,rng.get(),&pop.mut_lookup,generation,
						    0.005,0.,[&rng](){return gsl_rng_uniform(rng.get());},[](){return 0.;},[](){return 0.;}),
					  std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,std::placeholders::_4,
						    std::ref(pop.neutral),std::ref(pop.selected),
						    &pop.gametes,
						    0., //no rec
						    rng.get(),
						    recmap),
					  std::bind(KTfwd::insert_at_end<poptype::gamete_t,poptype::glist_t>,std::placeholders::_1,std::placeholders::_2),
					  []( poptype::dipvector_t::const_iterator dip ) { return KTfwd::multiplicative_diploid()(dip->first,dip->second,2.); },
					  std::bind(KTfwd::mutation_remover(),std::placeholders::_1,2*pop.N));
      KTfwd::update_mutations(&pop.mutations,&pop.fixations,&pop.fixation_times,&pop.mut_lookup,generation,2*pop.N);
    }

  //This is the test--we will add "ID" numbers to each diploid at the end of the sim
  //This has zero relevance to any model, and is just a test.
  unsigned i=0;
  std::for_each(pop.diploids.begin(),pop.diploids.end(),
		[&i]( diploid_t & dip ) { dip.i = i++; } );

  //attempt to copy
  KTfwd::serialize s;
  s(pop,mwriter(),diploid_writer());

  /*
    This is where this test differs from above.
    We get a string from s, create a new KTfwd serialize,
    and then attempt to deserialize from that.
  */
  std::string str(s.buffer.str());
  KTfwd::serialize s2;
  s2.buffer.str(str);
  
  poptype pop2(0);
  KTfwd::deserialize()(pop2,s2,mreader(),diploid_reader());
  BOOST_REQUIRE(pop.mutations.size() == pop2.mutations.size());
  //   unsigned NN=0;
  // for(const auto & g : pop.gametes)
  //   {
  //     if(g.n)++NN;
  //   }
  // BOOST_REQUIRE(NN == pop2.gametes.size());
  BOOST_REQUIRE(pop.gametes.size() == pop2.gametes.size());
  BOOST_REQUIRE(pop.diploids.size() == pop2.diploids.size());
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
  i=0;
  for( auto d1 = pop.diploids.begin(),d2 = pop2.diploids.begin() ; d1 != pop.diploids.end() ; ++d1,++d2,++i )
    {
      BOOST_CHECK(d1->i == d2->i);//This is the serialized "extra" data!!
      BOOST_CHECK(d2->i==i);
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

BOOST_AUTO_TEST_CASE( singlepop_serialized_copy_construct_test )
{
  using poptype = KTfwd::singlepop_serialized<mutation_t,mwriter,mreader,diploid_t,diploid_writer,diploid_reader>;
  poptype pop(1000);

  KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);
  
  //Evolve for 10 generations
  std::function<double(void)> recmap = std::bind(gsl_rng_uniform,rng.get());
  for( unsigned generation= 0 ; generation < 10 ; ++generation )
    {
      double wbar = KTfwd::sample_diploid(rng.get(),
					  &pop.gametes,
					  &pop.diploids,
					  &pop.mutations,
					  1000,
					  0.005,
					  std::bind(KTfwd::infsites(),std::placeholders::_1,std::placeholders::_2,rng.get(),&pop.mut_lookup,generation,
						    0.005,0.,[&rng](){return gsl_rng_uniform(rng.get());},[](){return 0.;},[](){return 0.;}),
					  std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,std::placeholders::_4,
						    std::ref(pop.neutral),std::ref(pop.selected),
						    &pop.gametes,
						    0., //no rec
						    rng.get(),
						    recmap),
					  std::bind(KTfwd::insert_at_end<poptype::gamete_t,poptype::glist_t>,std::placeholders::_1,std::placeholders::_2),
					  []( poptype::dipvector_t::const_iterator dip ) { return KTfwd::multiplicative_diploid()(dip->first,dip->second,2.); },
					  std::bind(KTfwd::mutation_remover(),std::placeholders::_1,2*pop.N));
      KTfwd::update_mutations(&pop.mutations,&pop.fixations,&pop.fixation_times,&pop.mut_lookup,generation,2*pop.N);
    }

  //This is the test--we will add "ID" numbers to each diploid at the end of the sim
  //This has zero relevance to any model, and is just a test.
  unsigned i=0;
  std::for_each(pop.diploids.begin(),pop.diploids.end(),
		[&i]( diploid_t & dip ) { dip.i = i++; } );

  //attempt to copy
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
  i=0;
  for( auto d1 = pop.diploids.begin(),d2 = pop2.diploids.begin() ; d1 != pop.diploids.end() ; ++d1,++d2,++i )
    {
      BOOST_CHECK(d1->i == d2->i);
      BOOST_CHECK(d2->i==i);
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
  using poptype = KTfwd::singlepop_serialized<mutation_t,mwriter,mreader,diploid_t,diploid_writer,diploid_reader>;
  poptype pop(1000);

  KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);
    
  //Evolve for 10 generations
  std::function<double(void)> recmap = std::bind(gsl_rng_uniform,rng.get());
  for( unsigned generation= 0 ; generation < 10 ; ++generation )
    {
      double wbar = KTfwd::sample_diploid(rng.get(),
					  &pop.gametes,
					  &pop.diploids,
					  &pop.mutations,
					  1000,
					  0.005,
					  std::bind(KTfwd::infsites(),std::placeholders::_1,std::placeholders::_2,rng.get(),&pop.mut_lookup,generation,
						    0.005,0.,[&rng](){return gsl_rng_uniform(rng.get());},[](){return 0.;},[](){return 0.;}),
					  std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,std::placeholders::_4,
						    std::ref(pop.neutral),std::ref(pop.selected),
						    &pop.gametes,
						    0., //no rec
						    rng.get(),
						    recmap),
					  std::bind(KTfwd::insert_at_end<poptype::gamete_t,poptype::glist_t>,std::placeholders::_1,std::placeholders::_2),
					  []( poptype::dipvector_t::const_iterator dip ) { return KTfwd::multiplicative_diploid()(dip->first,dip->second,2.); },
					  std::bind(KTfwd::mutation_remover(),std::placeholders::_1,2*pop.N));
      KTfwd::update_mutations(&pop.mutations,&pop.fixations,&pop.fixation_times,&pop.mut_lookup,generation,2*pop.N);
    }

  //This is the test--we will add "ID" numbers to each diploid at the end of the sim
  //This has zero relevance to any model, and is just a test.
  unsigned i=0;
  std::for_each(pop.diploids.begin(),pop.diploids.end(),
		[&i]( diploid_t & dip ) { dip.i = i++; } );

  //attempt to copy
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
  i=0;
  for( auto d1 = pop.diploids.begin(),d2 = pop2.diploids.begin() ; d1 != pop.diploids.end() ; ++d1,++d2,++i )
    {
      BOOST_CHECK(d1->i == d2->i);
      BOOST_CHECK(d2->i == i);
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
