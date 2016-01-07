/*! \include diploid_ind.cc
  Simulate a single, finite Wright-Fisher diploid population with mutation, recombination, and no selection.

  This is an individual-based implementation of the gamete-based program diploid.cc

  This program illustrates many features of fwdpp:
  1.  Custom mutation classes
  2.  Implementing a mutation model (infinitely-many sites)
  3.  Iterating a population through its life cycle
  4.  Outputting a sample in "ms" format
*/
#include <iostream>
#include <vector>
#include <Sequence/SimData.hpp>
#include <fwdpp/diploid.hh>
//Pull mutation model from fwdpp's "sugar" layer  (@ref md_md_sugar)
#include <fwdpp/sugar/infsites.hpp>
//typedef mutation_with_age mtype;
using mtype = KTfwd::popgenmut;
#define SINGLEPOP_SIM
#include <common_ind.hpp>

int main(int argc, char ** argv)
{
  if (argc != 8)
    {
      std::cerr << "Too few arguments\n"
		<< "Usage: diploid_ind N theta rho ngens samplesize nreps seed\n";
      exit(10);
    }
  int argument=1;
  const unsigned N = unsigned(atoi(argv[argument++]));           //Number of diploids
  const double theta = atof(argv[argument++]);         //4*n*mutation rate.  Note: mutation rate is per REGION, not SITE!!
  const double rho = atof(argv[argument++]);           //4*n*recombination rate.  Note: recombination rate is per REGION, not SITE!!
  const unsigned ngens = unsigned(atoi(argv[argument++]));       //Number of generations to simulate
  const unsigned samplesize1 = unsigned(atoi(argv[argument++])); //Sample size to draw from the population
  int nreps = atoi(argv[argument++]);                  //Number of replicates to simulate
  const unsigned seed = unsigned(atoi(argv[argument++]));        //Random number seed

  const double mu = theta/double(4*N);                 //per-gamete mutation rate

  /*
    littler r is the recombination rate per region per generation.

    For individual simulation (UNLIKE GAMETE-BASED SIMS!!!),
    r = rho/(4N)
  */
  const double littler = rho/double(4*N);

  //Write the command line to stderr
  std::copy(argv,argv+argc,std::ostream_iterator<char*>(std::cerr," "));
  std::cerr << '\n';

  //Initiate random number generation system (via fwdpp/sugar/GSLrng_t.hpp)
  GSLrng r(seed);

  unsigned twoN = 2*N;

  //recombination map is uniform[0,1)
  std::function<double(void)> recmap = std::bind(gsl_rng_uniform,r.get());

  while(nreps--)
    {
      /*
	 Initialize a population of N diploids via KTfwd::singlepop (fwdpp/sugar/singlepop.hpp)
	 Starting in fwdpp 0.3.3, the sugar struct constructor may take a second argument,
	 which reserves memory for the containers used to hold the intermediate results of
	 recombination events.  The default is 100, which is probably too small for large-4Nu
	 simulations, so we use the call to std::max to select something better.
      */
      singlepop_t pop(N);
      pop.mutations.reserve(size_t(std::ceil(std::log(2*N)*theta+0.667*theta)));
      unsigned generation;
      double wbar;

      for( generation = 0; generation < ngens; ++generation )
      	{
      	  //Iterate the population through 1 generation
      	  wbar = KTfwd::sample_diploid(r.get(),
      				       pop.gametes,  //non-const reference to gametes
      				       pop.diploids, //non-const reference to diploids
      				       pop.mutations, //non-const reference to mutations
				       pop.mcounts,
      				       N,     //current pop size, remains constant
      				       mu,    //mutation rate per gamete
      				       /*
      					 The mutation model (KTfwd::infsites) will be applied by
					 sample_diploid in order to add mutations to gametes each generation.
      				       */
				       std::bind(KTfwd::infsites(),std::placeholders::_1,std::placeholders::_2,r.get(),std::ref(pop.mut_lookup),generation,
						 mu,0.,[&r](){return gsl_rng_uniform(r.get());},[](){return 0.;},[](){return 0.;}),
				       //The recombination policy includes the uniform crossover rate
      				       std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,
						 std::placeholders::_3,std::placeholders::_4,
						 //Pass as reference
						 std::ref(pop.neutral),std::ref(pop.selected),
						 //Bummer--have to add reference wrapper now...
						 std::ref(pop.gametes),std::ref(pop.mutations),
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
      				       std::bind(KTfwd::emplace_back<singlepop_t::gamete_t,singlepop_t::gvec_t>,std::placeholders::_1,std::placeholders::_2),
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
      				       std::bind(KTfwd::multiplicative_diploid(),std::placeholders::_1,std::placeholders::_2,
						 std::placeholders::_3,2.),
      				       /*
      					 For each gamete still extant after sampling,
      					 remove the pointers to any mutations that have
      					 been fixed in the population.

      					 For more complex models such as quantitative
      					 traits evolving towards an optimum, one may not
      					 with to remove fixations.  In that case,
      					 replace the line below with
      					 std::bind(KTfwd::mutation_remover(),std::placeholders::_1,twoN));

      					 Under multiplicative fitness and Wright-Fisher sampling
      					 proportional to relative fitness, fixed mutations
      					 are just a constant applied to everyone's fitness, so we
      					 can remove them, making the simulation faster, etc.
      				       */
				       std::bind(KTfwd::mutation_remover(),std::placeholders::_1,twoN));
	  KTfwd::update_mutations(pop.mutations,pop.fixations,pop.fixation_times,pop.mut_lookup,pop.mcounts,generation,twoN);
	  assert(KTfwd::check_sum(pop.gametes,twoN));
	}
      Sequence::SimData sdata;
      
      std::vector<std::pair<double,std::string> > mslike = KTfwd::ms_sample(r.get(),pop.mutations,pop.gametes,pop.diploids,samplesize1,true);
      
      if(!mslike.empty())
	{
	  sdata.assign(mslike.begin(),mslike.end());
	  std::cout << sdata << '\n';
	}
      else
	{
	  std::cout << "//\nsegsites: 0\n";
	}
    }
  return 0;
}
