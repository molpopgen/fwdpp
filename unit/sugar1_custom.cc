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
#include <fwdpp/sugar/serialization.hpp>

using mutation_t = KTfwd::popgenmut;
using mwriter = KTfwd::mutation_writer;
using mreader = KTfwd::mutation_reader<mutation_t>;

//Custom diploid type.
struct diploid_t : public KTfwd::tags::custom_diploid_t
{
  using first_type = std::size_t;
  using second_type = std::size_t;
  first_type first;
  second_type second;
  unsigned i;
  diploid_t() : first(first_type()),second(second_type()),i(std::numeric_limits<unsigned>::max()) {}
  diploid_t(first_type g1, first_type g2) : first(g1),second(g2),i(std::numeric_limits<unsigned>::max()) {}
  bool operator==(const diploid_t & rhs) const
  {
    return this->first == rhs.first &&
    this->second == rhs.second &&
    this->i == rhs.i;
  }
};

struct diploid_writer
{
  using result_type = void;
  template<typename itr,typename streamtype>
  inline result_type operator()( itr i, streamtype & o ) const
  {
    o.write( reinterpret_cast<const char *>(&i.i),sizeof(unsigned) );
  }
};

struct diploid_reader
{
  using result_type = void;
  template<typename itr,typename streamtype>
  inline result_type operator()( itr i, streamtype & in ) const
  {
    in.read( reinterpret_cast<char *>(&i.i),sizeof(unsigned) );
  }
};

using singlepop_t = KTfwd::singlepop<mutation_t,diploid_t>;

void simulate( singlepop_t & pop )
{
  KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);

  //Evolve for 10 generations
  std::function<double(void)> recmap = std::bind(gsl_rng_uniform,rng.get());
  for( unsigned generation= 0 ; generation < 10 ; ++generation )
    {
      double wbar = KTfwd::sample_diploid(rng.get(),
					  pop.gametes,
					  pop.diploids,
					  pop.mutations,
					  pop.mcounts,
					  1000,
					  0.005,
					  std::bind(KTfwd::infsites(),std::placeholders::_1,std::placeholders::_2,rng.get(),std::ref(pop.mut_lookup),generation,
						    0.005,0.,[&rng](){return gsl_rng_uniform(rng.get());},[](){return 0.;},[](){return 0.;}),
					  std::bind(KTfwd::poisson_xover(),rng.get(),0.005,0.,1.,
						    std::placeholders::_1,std::placeholders::_2,std::placeholders::_3),
					  []( const diploid_t & dip, const singlepop_t::gcont_t & gametes,
					      const singlepop_t::mcont_t & mutations) { return KTfwd::multiplicative_diploid()(gametes[dip.first],gametes[dip.second],mutations,2.); },
					  pop.neutral,pop.selected);
      KTfwd::update_mutations(pop.mutations,pop.fixations,pop.fixation_times,pop.mut_lookup,pop.mcounts,generation,2*pop.N);
    }
}

BOOST_AUTO_TEST_CASE( singlepop_sugar_custom_test1 )
{
  singlepop_t pop(1000);
  simulate(pop);

  singlepop_t pop2(pop);

  BOOST_CHECK_EQUAL(pop==pop2,true);
}

BOOST_AUTO_TEST_CASE( singlepop_sugar_custom_test2 )
{
  singlepop_t pop(1000);
  simulate(pop);

  KTfwd::serialize s;
  s(pop,mwriter(),diploid_writer());
  singlepop_t pop2(0);
  KTfwd::deserialize()(pop2,s,mreader(),diploid_reader());
  BOOST_CHECK_EQUAL(pop==pop2,true);
}

BOOST_AUTO_TEST_CASE( singlepop_sugar_custom_test3 )
{
  singlepop_t pop(1000);
  simulate(pop);

  singlepop_t pop2(std::move(pop));
  //Should be false b/c move will leave
  //pop's containers in a wacky state
  BOOST_CHECK_EQUAL(pop==pop2,false);
}

BOOST_AUTO_TEST_CASE( singlepop_sugar_custom_test4 )
{
  singlepop_t pop(1000);
  simulate(pop);

  singlepop_t pop2=std::move(pop);
  //Should be false b/c move will leave
  //pop's containers in a wacky state
  BOOST_CHECK_EQUAL(pop==pop2,false);
}
