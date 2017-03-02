/*
In addition to the usual fwdpp depdencies, we need boost.python.

To compile:
g++ -fPIC -Wall -W -O3 -I.. -I. `python-config --includes` -std=c++11 -c social_evol.cc
g++ -std=c++11 -shared -o social_evol.so social_evol.o -lboost_python -lboost_system  -lpython -lgsl -lgslcblas

To run:
python test_boost_python.py
*/
#include <limits>
#include <algorithm>
#include <boost/python.hpp>
#include <boost/container/list.hpp>
#include <boost/container/vector.hpp>
#include <boost/pool/pool_alloc.hpp>
#include <boost/unordered_set.hpp>
#include <boost/functional/hash.hpp>

//Main fwdpp library header
#include <fwdpp/diploid.hh>
//Include the necessary "sugar" components
#include <fwdpp/sugar/singlepop.hpp>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/sugar/serialization.hpp>
#include <fwdpp/sugar/infsites.hpp>

//FWDPP-related stuff
using mtype = KTfwd::mutation;


//We need to define a custom diploid genotype for our model
struct diploid_t 
{
  using first_type = std::size_t;
  using second_type = std::size_t;
  first_type first;
  second_type second;
  unsigned i;
  //constructors, etc.
  diploid_t() : first(first_type()),second(second_type()),i(std::numeric_limits<unsigned>::max()) {}
  //"perfect forwarding" constructor does not work with iterator from boost containers...
  //diploid_t(first_type && g1, first_type && g2) : first(std::forward(g1)),second(std::forward(g2)),i(numeric_limits<unsigned>::max()) {}
  diploid_t(first_type g1, first_type g2) : first(g1),second(g2),i(std::numeric_limits<unsigned>::max()) {}
  //The following constructors SHOULD be generated automagically by your compiler, so you don't have to:
  //(no idea what, if any, performance effect this may have.  Worst case is prob. the move constructor doesn't get auto-generated...
  //diploid_t( const diploid_t & ) = default;
  //diploid_t( diploid_t && ) = default;
  //diploid_t & operator=(const diploid_t &) = default;
};

//Define our our population type via KTfwd::sugar 
using poptype = KTfwd::singlepop<mtype,diploid_t>;
/*
  We will use a gsl_rng_mt19937 as our RNG.
  This type is implicitly convertible to gsl_rng *,
  and auto-handles the gsl_rng_free steps, etc.
 */
using GSLrng = KTfwd::GSLrng_t<KTfwd::GSL_RNG_MT19937>;

 /*! \brief Continuous Snowdrift Game from Doebeli, Hauert, and Killingback (2004, Science, 306:859--862)
    as implemented by Wakano and Lehmann (2014, J Theor Biol, 351:83--95)
    \param phenotypes phenotypes in population
    \param fitnesses fitnesses in population
    \param b1 linear benefit term
    \param b2 quadratic benefit term
    \param c1 linear cost term
    \param c2 quadratic cost term
    \return Mean payoff (fitness) from pairwise interactions in snowdrift game
    \ingroup fitness
  */
struct snowdrift_diploid
{
  using result_type = double;
  inline result_type operator()( const poptype::diploid_t & dip,
				 const poptype::gcont_t &,
				 const poptype::mcont_t &,
				 const std::vector<double> & phenotypes, 
				 //Pass const refs to avoid copying!
				 const double & b1, const double & b2, const double & c1, const double & c2) const
  {
    unsigned N = phenotypes.size();
    double zself = phenotypes[dip.i];
    result_type fitness = 0;
    for (unsigned j = 0; j < N; ++j)	  
      {
	if (dip.i != j)
	  {
	    double zpair = zself+phenotypes[j];
	    // Payoff function from Fig 1
	    double a = b1*zpair + b2*zpair*zpair - c1*zself - c2*zself*zself;
	    fitness += 1 + std::max(a, 0.0);
	  }
      }
    return fitness/double(N-1);
  }
};

//END FWDPP-related stuff

//Evolve the population for some amount of time with mutation and recombination
poptype evolve( GSLrng & rng,
		const unsigned & N,
		const unsigned & generations,
		const double & mu,
		const double & mu_del,
		const double & recrate,
		const double & b1, const double & b2, const double & c1, const double & c2)
{
  poptype pop(N);
  pop.mutations.reserve(std::ceil(std::log(2*N)*(4.*double(N)*(mu+mu_del))+0.667*(4.*double(N)*(mu+mu_del))));
  std::function<double(void)> recmap = std::bind(gsl_rng_uniform,rng.get()); //uniform crossover map
  std::vector<double> phenotypes(N);
  for( unsigned generation = 0 ; generation < generations ; ++generation )
    {
      //Fill phenotypes
      unsigned i = 0;

      for(auto & dip : pop.diploids )
	{ 
	  dip.i = i; 
	  phenotypes[i++] = KTfwd::additive_diploid()(pop.gametes[dip.first],pop.gametes[dip.second],pop.mutations,2.); 
	}

      double wbar = KTfwd::sample_diploid(rng.get(),
					  pop.gametes,
					  pop.diploids,
					  pop.mutations,
					  pop.mcounts,
					  N,
					  mu+mu_del,
					  std::bind(KTfwd::infsites(),std::placeholders::_1,std::placeholders::_2,rng.get(),std::ref(pop.mut_lookup),
						    mu,mu_del,[&rng](){return gsl_rng_uniform(rng.get());},
						    [&rng](){return 0.1*(0.5-gsl_rng_uniform(rng.get()));},
						    [](){return 2.;}),
					  std::bind(KTfwd::poisson_xover(),rng.get(),recrate,0.,1.,
						    std::placeholders::_1,std::placeholders::_2,std::placeholders::_3),
					  std::bind(snowdrift_diploid(),std::placeholders::_1,std::placeholders::_2,
						    std::placeholders::_3,std::cref(phenotypes),b1,b2,c1,c2),
					  pop.neutral,pop.selected);

      KTfwd::update_mutations(pop.mutations,pop.fixations,pop.fixation_times,pop.mut_lookup,pop.mcounts,generation,2*pop.N);
    }
  return pop;
}

using namespace boost::python;

//Calculate the site-frequency spectrum for a sample
boost::python::list sfs(GSLrng & rng,const poptype & pop,const unsigned & nsam)
{
  std::map<double,unsigned> mutfreqs;
  unsigned twoN = 2*pop.N;

  for( unsigned i = 0 ; i < nsam ; ++i )
    {
      //pick a random chrom (w/replacement...)
      unsigned chrom = unsigned(gsl_ran_flat(rng.get(),0.,double(twoN)));
      //get reference to that chrom from the individual
      auto & gamete = (chrom%2==0.) ? pop.gametes[pop.diploids[chrom/2].first] : pop.gametes[pop.diploids[chrom/2].second];
      for( auto & m : gamete.mutations )
	{
	  auto pos = pop.mutations[m].pos;
	  auto pos_itr = mutfreqs.find(pos);
	  if( pos_itr == mutfreqs.end() )
	    {
	      mutfreqs.insert(std::make_pair(pos,1));
	    }
	  else
	    {
	      pos_itr->second++;
	    }
	}
    }
  //Now, fill in the SFS, omitting positions that are fixed in the sample
  std::vector<unsigned> __rv(nsam-1,0u);
  for( const auto & __x : mutfreqs )
    {
      if (__x.second < nsam) __rv[__x.second-1]++;
    }
  boost::python::list rv;
  for( const auto & __x : __rv ) rv.append(__x);
  return rv;
}

//Now, we can expose the stuff to python
BOOST_PYTHON_MODULE(social_evol)
{
  //Expose the type based on fwdpp's "sugar" layer
  class_<poptype>("poptype",init<unsigned>())
    .def("clear",&poptype::clear)
    ;
  //Expose the GSL wrapper
  class_<GSLrng>("GSLrng",init<unsigned>())
    ;
  //Expose the function to run the model
  def("evolve",evolve);
  //And one to get the sfs of a sample
  def("sfs",sfs);
}

