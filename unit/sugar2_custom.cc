/*! \file sugar2.cc
  \ingroup unit 
  \brief Testing KTfwd::metapop with custom diploid type
*/
#define BOOST_TEST_MODULE sugarTest2_custom
#define BOOST_TEST_DYN_LINK 

#include <config.h>
#include <algorithm>
#include <boost/test/unit_test.hpp>
#include <fwdpp/diploid.hh>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/sugar/metapop.hpp>
#include <fwdpp/sugar/infsites.hpp>

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

size_t migpop(const size_t & source_pop, gsl_rng * r, const double & mig_prob)
{
  if( gsl_rng_uniform(r) < mig_prob )
    {
      return ! source_pop;
    }
  return source_pop;
}

using poptype = KTfwd::metapop<mutation_t,diploid_t>;

void simulate(poptype & pop)
{
  //Evolve for 10 generations
  std::vector<std::function<double (const poptype::diploid_t &,
				    const poptype::gcont_t &,
				    const poptype::mcont_t &)> > fitness_funcs(2,
									       [](const poptype::diploid_t & d,
										  const poptype::gcont_t & g,
										  const poptype::mcont_t & m)
									       {
										 return KTfwd::multiplicative_diploid()(g[d.first],g[d.second],m);
									       }
									       );
  KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);
  for( unsigned generation= 0 ; generation < 10 ; ++generation )
    {
      std::vector<double> wbar = KTfwd::sample_diploid(rng.get(),
						       pop.gametes,
						       pop.diploids,
						       pop.mutations,
						       pop.mcounts,
						       &pop.Ns[0],
						       0.005,
						       std::bind(KTfwd::infsites(),std::placeholders::_1,std::placeholders::_2,rng.get(),std::ref(pop.mut_lookup),generation,
								 0.005,0.,[&rng](){return gsl_rng_uniform(rng.get());},[](){return 0.;},[](){return 0.;}),
						       std::bind(KTfwd::poisson_xover(),rng.get(),0.005,0.,1.,
								 std::placeholders::_1,std::placeholders::_2,std::placeholders::_3),
						       fitness_funcs,
						       std::bind(migpop,std::placeholders::_1,rng.get(),0.001),
						       pop.neutral,pop.selected);
      KTfwd::update_mutations(pop.mutations,pop.fixations,pop.fixation_times,pop.mut_lookup,pop.mcounts,generation,4000);
    }
}

BOOST_AUTO_TEST_CASE( metapop_sugar_custom_test1 )
{
  poptype pop({1000,1000});
  simulate(pop);
  poptype pop2(pop);
  BOOST_CHECK_EQUAL(pop==pop2,true);
}

BOOST_AUTO_TEST_CASE( metapop_sugar_custom_test2 )
{
  poptype pop({1000,1000});
  simulate(pop);
  poptype pop2{0,0};
  KTfwd::serialize s;
  s(pop,mwriter(),diploid_writer());
  KTfwd::deserialize()(pop2,s,mreader(),diploid_reader());
  BOOST_CHECK_EQUAL(pop==pop2,true);
}

BOOST_AUTO_TEST_CASE( metapop_sugar_custom_test3 )
{
  poptype pop({1000,1000});
  simulate(pop);
  poptype pop2(std::move(pop));
  BOOST_CHECK_EQUAL(pop==pop2,false);
}

BOOST_AUTO_TEST_CASE( metapop_sugar_custom_test4 )
{
  poptype pop({1000,1000});
  simulate(pop);
  poptype pop2=std::move(pop);
  BOOST_CHECK_EQUAL(pop==pop2,false);
}
