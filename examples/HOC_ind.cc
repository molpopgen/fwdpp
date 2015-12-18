/*! \include HOC_ind.cc
  A twist of additive House-of-Cards models with mutations within haplotypes
  having their effect sizes constrained to = new haplotype effect size.
*/
#include <config.h>
#include <vector>
#include <list>
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
  inline singlepop_t::mlist_t::iterator
  operator()( mqueue_t & mut_recycling_bin,
	      poptype::gamete_t & g,
	      singlepop_t::mlist_t * mutations,
	      gsl_rng * r, poptype::lookup_table_t * lookup, const double & sigmu) const
  {
    double pos = gsl_rng_uniform(r);
    while( lookup->find(pos) != lookup->end() ) //make sure it doesn't exist in the population
      { 
    	pos = gsl_rng_uniform(r);  //if it does, generate a new one
      }
    lookup->insert(pos);
    double E = gsl_ran_gaussian(r,sigmu); //effect size of hap, after mutation
    double sum = std::accumulate(g.smutations.begin(),g.smutations.end(),0.,
     				 [](double & d, const poptype::mlist_t::iterator & m) { return d + m->s; });
    double esize = (E > sum) ? std::fabs(E-sum) : -1.*std::fabs(E-sum);
    if(!mut_recycling_bin.empty())
      {
	auto rv = mut_recycling_bin.front();
	mut_recycling_bin.pop();
	rv->pos=pos;
	rv->s=esize;
	rv->n=1u;
	return rv;
      }
    return mutations->emplace(mutations->end(),pos,esize,1);
    // return mtype(pos,esize,1);
  }
};

double simple_gaussian( const poptype::glist_t::iterator & h1, const poptype::glist_t::iterator & h2 )
{
  double sum = std::accumulate(h1->smutations.begin(),h1->smutations.end(),0.,
			       [](double & d, const poptype::mlist_t::iterator & m) { return d + m->s; });
  sum += std::accumulate(h2->smutations.begin(),h2->smutations.end(),0.,
			 [](double & d, const poptype::mlist_t::iterator & m) { return d + m->s; });
  return std::exp(-1.*std::pow(sum,2.)/2.);
}

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
      exit(10);
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

  unsigned twoN = 2*N;

  //recombination map is uniform[0,1)
  std::function<double(void)> recmap = std::bind(gsl_rng_uniform,r.get());

  while(nreps--)
    {
      poptype pop(N);

      //HOChap mmodel(r,sigmu,&pop.mut_lookup);
      for( unsigned generation = 0; generation < ngens; ++generation )
      	{
      	  //Iterate the population through 1 generation
      	  double wbar = KTfwd::sample_diploid(r.get(),
					      &pop.gametes,  //non-const pointer to gametes
					      &pop.diploids, //non-const pointer to diploids
					      &pop.mutations, //non-const pointer to mutations
					      N,     //current pop size, remains constant
					      mu,    //mutation rate per gamete
					      /*
						The mutation model below will be applied by
						sample_diploid in order to add mutations to gametes each generation.

						The _1 is a placeholder for a non-const reference to a gamete (see defn'
						of HOChap above).
					      */
					      std::bind(HOChap(),std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,
							r.get(),&pop.mut_lookup,sigmu),
					      // mmodel,
					      //The recombination policy includes the uniform crossover rate
					      std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,std::placeholders::_4,
							std::ref(pop.neutral),std::ref(pop.selected),
							&pop.gametes,
							littler,
							r.get(),
							recmap),
					      /*
						Policy telling KTfwd::mutate how to add mutated gametes into the gamete pool.
						If mutation results in a new gamete, add that gamete to the 
						end of gametes. This is always the case under infinitely-many sites,
						but for other mutation models, mutation may result in a new
						copy identical to an existing gamete.  If so,
						that gamete's frequency increases by 1.
					      */
					      std::bind(KTfwd::insert_at_end<poptype::gamete_t,poptype::glist_t>,std::placeholders::_1,std::placeholders::_2),
					      /*
						Fitness is multiplicative over sites.
						
						The fitness model takes two gametes as arguments.  
						The gametes are passed to this function by 
						KTfwd::sample_diploid, and the _1 and _2 are placeholders for
						those gametes (see documentation for boost/bind.hpp for details).
						The 2. means that fitnesses are 1, 1+sh, and 1+2s for genotypes
						AA, Aa, and aa, respectively, where a is a mutation with
						selection coefficient s and dominance h, and the fitness of 
						the diploid is the product of fitness over sites
						
						There is no selection in this simulation, but this
						function is called anyways to illustrate it as multiplicative
						models are very common in population genetics
					      */
					      std::bind(simple_gaussian,std::placeholders::_1,std::placeholders::_2),
					      KTfwd::remove_nothing());
      	  KTfwd::update_mutations(&pop.mutations,&pop.mut_lookup);
	  assert(KTfwd::check_sum(pop.gametes,twoN));
	}    
      //Get VG for this replicate
      //Note: this can be done more efficiently via the boost accumulator library, which
      //we don't use here to reduce dependencies.
      std::vector<double> G;
      std::for_each(pop.diploids.cbegin(),pop.diploids.cend(),[&G](const poptype::diploid_t & dip ) {
      	  double sum = std::accumulate(dip.first->smutations.begin(),dip.first->smutations.end(),0.,
      				       [](double & d, const poptype::mlist_t::iterator & m) { return d + m->s; } );
      	  sum += std::accumulate(dip.second->smutations.begin(),dip.second->smutations.end(),0.,
      				 [](double & d, const poptype::mlist_t::iterator & m) { return d + m->s; } );
	  G.push_back(sum);
	});
      std::cout << gsl_stats_variance(&G[0],1,G.size()) << std::endl;
    }
  return 0;
}

