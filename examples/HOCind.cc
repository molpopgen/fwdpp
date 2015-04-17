/*! \include HOCind.cc
  A twist of additive House-of-Cards models with mutations within haplotypes
  having their effect sizes constrained to = new haplotype effect size.
*/

#include <fwdpp/diploid.hh>
#include <Sequence/SimData.hpp>
#include <vector>
#include <list>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/variance.hpp>
/*
  NOTE: we are using boost containers with boost memory allocators.  
  My experience on my systems is that each of those objects is worth a 5-10% 
  run time speedup compare to the same objects from namespace std.

  fwdpp is compatible with ANY container system as long as the
  one containing gametes conforms to the behavior of std::vector and 
  the one containing mutations conforms to the behavior of std::list 
  (especially w/regards to no pointer invalidation upon insertion/delete!!)
*/

//compiling the code with -DUSE_STANDARD_CONTAINERS will use std::vector and std::list instead of the boost alternatives
typedef KTfwd::mutation mtype;
#include <common_ind.hpp>

struct HOChap : public KTfwd::tags::gamete_dependent
{
  gsl_rng * r;
  const double sigmu;
  lookup_table_type * lookup;
  HOChap( gsl_rng * __r, const double & __sigmu,
	  lookup_table_type * __lookup ) : r(__r),sigmu(__sigmu),lookup(__lookup)
  {
  }
  inline mtype operator()(gtype & g,mlist * mutations) const
  {
    double pos = gsl_rng_uniform(r);
    while( lookup->find(pos) != lookup->end() ) //make sure it doesn't exist in the population
      { 
	pos = gsl_rng_uniform(r);  //if it does, generate a new one
      }
    lookup->insert(pos);
    double E = gsl_ran_gaussian(r,sigmu); //effect size of hap, after mutation
    double sum = std::accumulate(g.smutations.begin(),g.smutations.end(),0.,
				 [](double & d, const mlist::iterator & m) { return d + m->s; });
    double esize = (E > sum) ? std::fabs(E-sum) : -1.*std::fabs(E-sum);
    return mtype(pos,esize,1);
  }
};

double simple_gaussian( const glist::iterator & h1, const glist::iterator & h2 )
{
  double sum = std::accumulate(h1->smutations.begin(),h1->smutations.end(),0.,
			       [](double & d, const mlist::iterator & m) { return d + m->s; });
  sum += std::accumulate(h2->smutations.begin(),h2->smutations.end(),0.,
			 [](double & d, const mlist::iterator & m) { return d + m->s; });
  return std::exp(-1.*std::pow(sum,2.)/2.);
}

int main(int argc, char ** argv)
{
  if (argc != 8)
    {
      std::cerr << "Too few arguments\n"
		<< "Usage: diploid_ind N theta rho ngens samplesize nreps seed\n";
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
  
  //Initiate random number generation system
  gsl_rng * r =  gsl_rng_alloc(gsl_rng_ranlxs2);
  gsl_rng_set(r,seed);

  unsigned twoN = 2*N;

  //recombination map is uniform[0,1)
  std::function<double(void)> recmap = std::bind(gsl_rng_uniform,r);

  while(nreps--)
    {
      //the population begins with 1 gamete with no mutations amd initial count 2N
      glist gametes(1,gtype(twoN));

      //Establish the list of diploids.  Here, there is only 1 gamete, so it is easy:
      std::vector< std::pair< glist::iterator, glist::iterator> > diploids(N,
									   std::make_pair(gametes.begin(),
											  gametes.begin()));
      mlist mutations;  //the population is devoid of mutations
      std::vector<mtype> fixations;  //store mutation that fix.  Passed to KTfwd::remove_fixed_lost
      std::vector<unsigned> fixation_times; //store when they fix.  Passed to KTfwd::remove_fixed_lost
      unsigned generation;
      double wbar;
      lookup_table_type lookup;  //this is our lookup table for the mutation model

      HOChap mmodel(r,sigmu,&lookup);
      for( generation = 0; generation < ngens; ++generation )
      	{
      	  //Iterate the population through 1 generation
      	  wbar = KTfwd::sample_diploid(r,
      				       &gametes,  //non-const pointer to gametes
      				       &diploids, //non-const pointer to diploids
      				       &mutations, //non-const pointer to mutations
      				       N,     //current pop size, remains constant
      				       mu,    //mutation rate per gamete
      				       /*
      					 The mutation model (defined above) will pass each gamete
      					 to be mutated to the mutation model function.  Again, _1
      					 is used as a placeholder for that gamete.
      				       */
				       mmodel,
				       //The recombination policy includes the uniform crossover rate
      				       std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,
						   &gametes,
      						   littler,
      						   r,
      						   recmap),
				       /*
					 Policy to insert new mutations at the end of the mutations list
				       */
      				       std::bind(KTfwd::insert_at_end<mtype,mlist>,std::placeholders::_1,std::placeholders::_2),
				       /*
					 Policy telling KTfwd::mutate how to add mutated gametes into the gamete pool.
					 If mutation results in a new gamete, add that gamete to the 
					 end of gametes. This is always the case under infinitely-many sites,
					 but for other mutation models, mutation may result in a new
					 copy identical to an existing gamete.  If so,
					 that gamete's frequency increases by 1.
				       */
      				       std::bind(KTfwd::insert_at_end<gtype,glist>,std::placeholders::_1,std::placeholders::_2),
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
      				       /*
      					 For each gamete still extant after sampling,
      					 remove the pointers to any mutations that have 
      					 been lost from the population.
      				       */
      				       std::bind(KTfwd::mutation_remover(),std::placeholders::_1,0));
      	  KTfwd::remove_lost(&mutations,&lookup);
	  assert(KTfwd::check_sum(gametes,twoN));
	}    
      //Get VG for this replicate
      boost::accumulators::accumulator_set<double,boost::accumulators::stats<boost::accumulators::tag::variance> > VG;
      std::for_each(diploids.cbegin(),diploids.cend(),[&VG](const std::pair< glist::iterator, glist::iterator> & dip ) {
	  double sum = std::accumulate(dip.first->smutations.begin(),dip.first->smutations.end(),0.,
				       [](double & d, const mlist::iterator & m) { return d + m->s; } );
	  sum += std::accumulate(dip.second->smutations.begin(),dip.second->smutations.end(),0.,
				 [](double & d, const mlist::iterator & m) { return d + m->s; } );
	  VG(sum);
	} );
      std::cout << boost::accumulators::variance(VG) << '\n';
    }
  gsl_rng_free(r);
  return 0;
}
