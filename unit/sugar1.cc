//! \file sugar1.cc \ingroup unit
#define BOOST_TEST_MODULE sugarTest1
#define BOOST_TEST_DYN_LINK 

#include <config.h>
#include <boost/test/unit_test.hpp>
#include <gsl/gsl_rng.h>
#include <fwdpp/diploid.hh>
#define FWDPP_SUGAR_USE_BOOST
#include <fwdpp/sugar/singlepop.hpp>

struct mutation_with_age : public KTfwd::mutation_base
{
  unsigned g;
  double s,h;
  /*
    The constructor initializes all class data, including that of the base class via a constructor
    call to the base class.
  */
  mutation_with_age(const unsigned & __o,const double & position, const unsigned & count, const bool & isneutral = true)
    : KTfwd::mutation_base(position,count,isneutral),g(__o),s(0.),h(0.)
  {	
  }
};

//function object to write mutation data in binary format
struct mwriter
{
  typedef void result_type;
  result_type operator()( const mutation_with_age & m, std::ostringstream & buffer ) const
  {
    unsigned u = m.n;
    buffer.write( reinterpret_cast< char * >(&u),sizeof(unsigned) );
    u = m.g;
    buffer.write( reinterpret_cast< char * >(&u),sizeof(unsigned) );
    bool b = m.neutral;
    buffer.write( reinterpret_cast< char * >(&b),sizeof(bool) );
    double d = m.pos;
    buffer.write( reinterpret_cast< char * >(&d),sizeof(double) );
    d = m.s;
    buffer.write( reinterpret_cast< char * >(&d),sizeof(double) );
  }
};

//function object to read mutation data in binary format
struct mreader
{
  typedef mutation_with_age result_type;
  result_type operator()( std::istream & in ) const
  {
    unsigned n;
    in.read( reinterpret_cast< char * >(&n),sizeof(unsigned) );
    unsigned g;
    in.read( reinterpret_cast< char * >(&g),sizeof(unsigned) );
    bool neut;
    in.read( reinterpret_cast< char * >(&neut),sizeof(bool) );
    double pos;
    in.read( reinterpret_cast< char * >(&pos),sizeof(double) );
    double s;
    in.read( reinterpret_cast< char * >(&s),sizeof(double) );
    return result_type(g,pos,n);
  }
};

using poptype = KTfwd::singlepop<mutation_with_age>;

mutation_with_age neutral_mutations_inf_sites(gsl_rng * r,const unsigned & generation,poptype::mlist_t * mutations,
					      KTfwd::singlepop<mutation_with_age>::lookup_table_t * lookup)
{
  double pos = gsl_rng_uniform(r);
  while( lookup->find(pos) != lookup->end() ) //make sure it doesn't exist in the population
    { 
      pos = gsl_rng_uniform(r);  //if it does, generate a new one
    }
  lookup->insert(pos);
  return mutation_with_age(generation,pos,1,true);
}

BOOST_AUTO_TEST_CASE( test1 )
{
  poptype pop(1000);

  gsl_rng * r =  gsl_rng_alloc(gsl_rng_ranlxs2);
  gsl_rng_set(r,0); //seed is set to 0

  //Evolve for 10 generations
  std::function<double(void)> recmap = std::bind(gsl_rng_uniform,r);
  for( unsigned generation= 0 ; generation < 10 ; ++generation )
    {
      double wbar = KTfwd::sample_diploid(r,
					  &pop.gametes,
					  &pop.diploids,
					  &pop.mutations,
					  1000,
					  0.005,
					  std::bind(neutral_mutations_inf_sites,r,generation,std::placeholders::_1,&pop.mut_lookup),
					  std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,
						    &pop.gametes,
						    0., //no rec
						    r,
						    recmap),
					  std::bind(KTfwd::insert_at_end<poptype::mutation_t,poptype::mlist_t>,std::placeholders::_1,std::placeholders::_2),
					  std::bind(KTfwd::insert_at_end<poptype::gamete_t,poptype::glist_t>,std::placeholders::_1,std::placeholders::_2),
					  std::bind(KTfwd::multiplicative_diploid(),std::placeholders::_1,std::placeholders::_2,2.),
					  std::bind(KTfwd::mutation_remover(),std::placeholders::_1,0,2*pop.N));
      KTfwd::remove_fixed_lost(&pop.mutations,&pop.fixations,&pop.fixation_times,&pop.mut_lookup,generation,2*pop.N);
    }

  //attempt to copy
  KTfwd::serialize s;
  s(pop,mwriter());
  poptype pop2(0);
  KTfwd::deserialize()(pop2,s,mreader());

  BOOST_REQUIRE(pop.mutations.size() == pop2.mutations.size());
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
