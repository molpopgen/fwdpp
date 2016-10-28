/*! \include HOC_ind.cc
  A twist of additive House-of-Cards models with mutations within haplotypes
  having their effect sizes constrained to = new haplotype effect size.
*/
#include <config.h>
#include <vector>
#include <iostream>
#include <gsl/gsl_statistics_double.h>
#include <fwdpp/diploid.hh>
#include <fwdpp/sugar/GSLrng_t.hpp>

using mtype = KTfwd::mutation;
#define SINGLEPOP_SIM
#include <common_ind.hpp>
using poptype = singlepop_t;

struct HOChap
{
  template<typename mqueue_t>
  std::size_t
  operator()( mqueue_t & mut_recycling_bin,
	      poptype::gamete_t & g,
	      singlepop_t::mcont_t & mutations,
	      const gsl_rng * r, poptype::lookup_table_t & lookup, const double & sigmu) const
  {
    double pos = gsl_rng_uniform(r);
    while( lookup.find(pos) != lookup.end() ) //make sure it doesn't exist in the population
      { 
    	pos = gsl_rng_uniform(r);  //if it does, generate a new one
      }
    lookup.insert(pos);
    double E = gsl_ran_gaussian(r,sigmu); //effect size of hap, after mutation
    double sum = std::accumulate(g.smutations.begin(),g.smutations.end(),0.,
     				 [&mutations](double & d, const std::size_t m) { return d + mutations[m].s; });
    double esize = (E > sum) ? std::fabs(E-sum) : -1.*std::fabs(E-sum);
    return KTfwd::fwdpp_internal::recycle_mutation_helper(mut_recycling_bin,mutations,pos,esize,1);
  }
};

double addEsizes(const poptype::gamete_t & h,
		 const poptype::mcont_t & mutations)
{
  return std::accumulate(h.smutations.begin(),h.smutations.end(),0.,
			 [&mutations](double & d, const std::size_t & m) { return d + mutations[m].s; });
}

double gaussianFitness(const double & a, const double & b) {return std::exp(-1.*std::pow(a+b,2.)/2.); }

int main(int argc, char ** argv)
{
  if (argc != 8)
    {
      std::cerr << "Too few arguments\n"
		<< "Usage: " << argv[0] << " N mu sigmu rho ngens nreps seed\n"
		<< "where:\n"
		<< "N = diploid population size\n"
		<< "mu = mutation rate to variants affecting fitness\n"
		<< "sigmu = standard deviation of effect sizes.  E ~ N(0,sigmu^2)\n"
		<< "4Nr = population scaledrecombination rate (per diploid, per generation)\n"
		<< "ngens = number of generations to simulate\n"
		<< "nreps = number of replicates to simulate\n"
		<< "seed = random number seed\n";	
      exit(0);
    } 
  int argument=1;
  const unsigned N = atoi(argv[argument++]);           //Number of diploids
  const double mu = atof(argv[argument++]);            //mutation rate to mutations affecting trait
  const double sigmu = atof(argv[argument++]);         //std. dev. of effect size
  const double rho = atof(argv[argument++]);           //4*n*recombination rate.  Note: recombination rate is per REGION, not SITE!!
  const unsigned ngens = atoi(argv[argument++]);       //Number of generations to simulate
  int nreps = atoi(argv[argument++]);                  //Number of replicates to simulate
  const unsigned seed = atoi(argv[argument++]);        //Random number seed

  /*
    littler r is the recombination rate per region per generation.

    For individual simulation (UNLIKE GAMETE-BASED SIMS!!!),
    r = rho/(4N)
  */
  const double littler = rho/double(4*N);
  
  //Write the command line to stderr
  std::copy(argv,argv+argc,std::ostream_iterator<char*>(std::cerr," "));
  std::cerr << '\n';
  
  //Initiate random number generation system from sugar layer
  KTfwd::GSLrng_t<KTfwd::GSL_RNG_MT19937> r(seed);

  //recombination map is uniform[0,1)
  std::function<double(void)> recmap = std::bind(gsl_rng_uniform,r.get());

  while(nreps--)
    {
      poptype pop(N);
      pop.mutations.reserve(std::ceil(std::log(2*N)*(4.*double(N)*mu)+0.667*(4.*double(N)*mu)));
      //HOChap mmodel(r,sigmu,&pop.mut_lookup);
      for( unsigned generation = 0; generation < ngens; ++generation )
      	{
      	  //Iterate the population through 1 generation
      	  double wbar = KTfwd::sample_diploid(r.get(),
					      pop.gametes,  //non-const pointer to gametes
					      pop.diploids, //non-const pointer to diploids
					      pop.mutations, //non-const pointer to mutations
					      pop.mcounts,
					      N,     //current pop size, remains constant
					      mu,    //mutation rate per gamete
					      /*
						The mutation model below will be applied by
						sample_diploid in order to add mutations to gametes each generation.

						The _1 is a placeholder for a non-const reference to a gamete (see defn'
						of HOChap above).
					      */
					      std::bind(HOChap(),std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,
							r.get(),std::ref(pop.mut_lookup),sigmu),
					      std::bind(KTfwd::poisson_xover(),r.get(),littler,0.,1.,
						 std::placeholders::_1,std::placeholders::_2,std::placeholders::_3),
					      std::bind(KTfwd::haplotype_dependent_fitness(),
							std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,
							addEsizes,gaussianFitness),
					      pop.neutral,pop.selected,
					      0.,
					      KTfwd::remove_nothing());
      	  KTfwd::update_mutations(pop.mutations,pop.mut_lookup,pop.mcounts,2*N);
	  assert(KTfwd::check_sum(pop.gametes,2*N));
	}    
      //Get VG for this replicate
      //Note: this can be done more efficiently via the boost accumulator library, which
      //we don't use here to reduce dependencies.
      std::vector<double> G;
      for(const auto & d : pop.diploids)
	{
	  double sum = 0.;
	  for(const auto & m : pop.gametes[d.first].smutations) sum += pop.mutations[m].s;
	  for(const auto & m : pop.gametes[d.second].smutations) sum += pop.mutations[m].s;
	  G.push_back(sum);
	}
      std::cout << gsl_stats_variance(&G[0],1,G.size()) << std::endl;
    }
  return 0;
}

